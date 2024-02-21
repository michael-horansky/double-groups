# ---------------------------------------------------------
#
#   The leading src file of the AReTDoG package
#   AReTDoG: Algorithmic Representation Theory of Double Groups
#   Import this file to access the full scope of AReTDoG
#
#   This file provides:
#       class QDGroup(Group)
# ---------------------------------------------------------
#   Created: 8 November 2023
#   Author: Michal Horansky (michal.horansky20@imperial.ac.uk)
#
# ---------------------------------------------------------


import tikz
import os
from .class_group import *

# hardcoded solutions yay!!!
tex_readable_irrep_labels = {"D3h QD" : {
        "A_1" : "A_1'",
        "A_2" : "A_2''",
        "A_3" : "A_1''",
        "A_4" : "A_2'",
        "E_1" : "E''",
        "E_2" : "E'",
        "E_1(j=3/2)" : "E_{3/2}",
        "E_2(j=1/2)" : "E_{1/2}",
        "E_3(j=5/2)" : "E_{5/2}"
    }, "C3v QD" : {
        "A_1" : "A_1",
        "A_2" : "A_2",
        "E_1" : "E",
        "E_1(j=1/2)" : "E_{1/2}",
        "1E(j=3/2)" : "{}^1E_{3/2}",
        "2E(j=3/2)" : "{}^2E_{3/2}"
    }}

tikz_line_thickness = ["very thin", "thick", "ultra thick"]


class QDGroup(Group):
    # Instances of this class possess properties specific to the semiconductor structure, like the conduction and the valence band orbitals, the notion of holes etc
    
    def __init__(self, name):
        self.holes = {} # {"hole name" : Rep}
        self.hole_irreps = {} # {"hole name" : "hole irrep or duo of conjugate irreps"}
        self.electrons = {} # we can have multiple characters of electrons when in higher orbital (very excited)
        self.electron_irreps = {} # "rep name"
        self.excitons = {} # {"label" : Rep}
        self.exciton_irreps = {} # {"label" : "rep name"}
        self.exciton_occupancies = {} # {"label" : [[electron occupancies], [hole_occupancies]]}
        
        self.transition_chain = {} # {"exciton label" : [ ["exciton label", recombining_e_index, recombining_h_index], ... ]}
        
        super().__init__(name)
    
    
    
    # --------------- exciton labelling methods
    
    
    def classify_fermions(self, orbital, spin=1/2, band="valence"):
        # for GaAs: conduction band usually has s-orbital and valence band usually has p-orbital
        # if band="valence", we select highest possible j (highest energy). If "conduction", we select lowest j
        
        orbital_dict = {"s" : 0, "p" : 1, "d" : 2, "f" : 3}
        if type(orbital) == str: #otherwise must be a number
            orbital = orbital_dict[orbital]
        if band == "valence" or band == "v":
            j_val = orbital + spin
        elif band == "conduction" or band == "c":
            j_val = np.abs(orbital - spin)
        gerade_rep = self.angular_representation(j_val, "g")
        if orbital % 2 == 1:
            gerade_rep *= self.irrep_characters[self.inversion_irrep]
        
        # now we classify these:
        return(self.separate_constituent_representations(self.reduce_representation(gerade_rep)[0]))
    
    def classify_electrons(self, orbital):
        electron_dict = self.classify_fermions(orbital, 1/2, "conduction")
        electron_i = 1
        for irrep_label, rep in electron_dict.items():
            self.electrons[f"e{electron_i}"] = rep
            self.electron_irreps[f"e{electron_i}"] = self.reduce_representation(rep)[1]
            electron_i += 1
    
    def classify_holes(self, orbital):
        hole_dict = self.classify_fermions(orbital, 1/2, "valence")
        hole_i = 1
        for irrep_label, rep in hole_dict.items():
            self.holes[f"h{hole_i}"] = rep
            self.hole_irreps[f"h{hole_i}"] = self.reduce_representation(rep)[1]
            hole_i += 1
    
    
    def exciton_rep(self, electron_numbers, hole_numbers):
        result_rep = self.irrep_characters[self.identity_irrep] #identity rep
        i = 0
        for e_name, e_rep in self.electrons.items():
            # pauli exclusion - each filled level is labelled by identity
            if type(electron_numbers) == dict:
                free_e = electron_numbers[e_name] % int(np.real(e_rep.characters["E"]))
            else:
                free_e = electron_numbers[i] % int(np.real(e_rep.characters["E"]))
            
            # what to do when multiple free e? idk mate
            for j in range(free_e):
                result_rep *= e_rep
            i += 1
        i = 0
        for h_name, h_rep in self.holes.items():
            # pauli exclusion - each filled level is labelled by identity
            if type(hole_numbers) == dict:
                free_h = hole_numbers[h_name] % int(np.real(h_rep.characters["E"]))
            else:
                free_h = hole_numbers[i] % int(np.real(h_rep.characters["E"]))
            
            # what to do when multiple free e? idk mate
            for j in range(free_h):
                result_rep *= h_rep
            i += 1
        
        return(result_rep)
    
    
    def classify_excitons(self, max_e, max_h, e_orbital = "s", h_orbital = "p"):
        # max_e is the maximum amount of total electrons
        # max_h is the maximum amount of total holes
        
        # Since we allow multiple electron characters to exist, we generalise the notation to (number of e-h pairs)X[e1, e2, e3...][h1, h2, h3...]
        # if number of e-h pairs = 0: the notation becomes a h_a + b h_b +...; a e_a + b e_b +..., or vacuum for the trivial case
        
        self.classify_electrons(e_orbital)
        self.classify_holes(h_orbital)
        
        electron_labels = list(self.electrons.keys())
        hole_labels = list(self.holes.keys())
        
        N_e = len(electron_labels)
        N_h = len(hole_labels)
        
        # first, vacuum
        self.excitons = {"vacuum" : self.irrep_characters[self.identity_irrep]}
        self.exciton_irreps = {"vacuum" : self.identity_irrep}
        self.exciton_occupancies = {"vacuum" : [[0] * N_e, [0] * N_h]}
        
        
        e_occ_by_e_n = [[[0] * N_e]]
        h_occ_by_h_n = [[[0] * N_h]]
        
        
        # now, for hole = 0, e != 0:
        hole_occupancies = [0] * N_h
        for electron_number in range(1, max_e + 1):
            electron_occupancies = integer_partitions(electron_number, N_e)
            e_occ_by_e_n.append(electron_occupancies)
            for occupancy in electron_occupancies:
                label = ""
                for i in range(N_e):
                    if occupancy[i] > 0:
                        label += f"{occupancy[i]}{electron_labels[i]}+"
                label = label[:-1]
                self.excitons[label] = self.exciton_rep(occupancy, hole_occupancies)
                self.exciton_irreps[label] = self.reduce_representation(self.excitons[label])[1]
                self.exciton_occupancies[label] = [occupancy, hole_occupancies]
        
        # now, for electron = 0, hole != 0:
        electron_occupancies = [0] * N_e
        for hole_number in range(1, max_h + 1):
            hole_occupancies = integer_partitions(hole_number, N_h)
            h_occ_by_h_n.append(hole_occupancies)
            for occupancy in hole_occupancies:
                label = ""
                for i in range(N_h):
                    if occupancy[i] > 0:
                        label += f"{occupancy[i]}{hole_labels[i]}+"
                label = label[:-1]
                self.excitons[label] = self.exciton_rep(electron_occupancies, occupancy)
                self.exciton_irreps[label] = self.reduce_representation(self.excitons[label])[1]
                self.exciton_occupancies[label] = [electron_occupancies, occupancy]
        
        # now, for the rest
        for electron_number in range(1, max_e + 1):
            for hole_number in range(1, max_h + 1):
                exciton_prefix = f"{min(electron_number, hole_number)}"
                for electron_occupancy in e_occ_by_e_n[electron_number]:
                    for hole_occupancy in h_occ_by_h_n[hole_number]:
                        label = exciton_prefix + "X[" + ",".join([str(x) for x in electron_occupancy]) + "][" + ",".join([str(x) for x in hole_occupancy]) + "]"
                        self.excitons[label] = self.exciton_rep(electron_occupancy, hole_occupancy)
                        self.exciton_irreps[label] = self.reduce_representation(self.excitons[label])[1]
                        self.exciton_occupancies[label] = [electron_occupancy, hole_occupancy]
                        
    def print_exciton_complexes(self):
        max_label_len = max([len(l) for l in self.exciton_irreps.keys()])
        for label, hr_rep in self.exciton_irreps.items():
            print(f"{' ' * (max_label_len - len(label)) + label} - {hr_rep}")
    
    
    def find_transition_chain(self):
        # by convention - if no transition from exciton, it is still included as a key, but its value is an empty list
        self.transition_chain = {}
        
        def recombined_exciton_complex(electron_occupancy, hole_occupancy, e_i, h_i):
            new_e_occ = []
            new_h_occ = []
            for i in range(len(electron_occupancy)):
                if i == e_i:
                    new_e_occ.append(electron_occupancy[i] - 1)
                else:
                    new_e_occ.append(electron_occupancy[i])
            for i in range(len(hole_occupancy)):
                if i == h_i:
                    new_h_occ.append(hole_occupancy[i] - 1)
                else:
                    new_h_occ.append(hole_occupancy[i])
            return([new_e_occ, new_h_occ])
                    
        
        
        for exciton_label in self.excitons.keys():
            # checks all non-zero occupancies for recombination
            self.transition_chain[exciton_label] = []
            cur_occupancies = self.exciton_occupancies[exciton_label]
            for i_e in range(len(cur_occupancies[0])):
                if cur_occupancies[0][i_e] > 0:
                    for i_h in range(len(cur_occupancies[1])):
                        if cur_occupancies[1][i_h] > 0:
                            # can recombine
                            new_occupancies = recombined_exciton_complex(cur_occupancies[0], cur_occupancies[1], i_e, i_h)
                            if new_occupancies in self.exciton_occupancies.values():
                                target_label = list(self.exciton_occupancies.keys())[list(self.exciton_occupancies.values()).index(new_occupancies)]
                                self.transition_chain[exciton_label].append(target_label)
                            else:
                                print(f"ERROR! Encountered a possible exciton recombination ({exciton_label}, {i_e}, {i_h}) that isn't classified under excitons!")
    
    # ---------------------------------------------------
    # --------------- OOUTPUT METHODS -------------------
    # ---------------------------------------------------
    
    def tex_readable_exciton_labels(self, x_label):
        #return(x_label)
        if x_label == "vacuum":
            return("vac.")
        # check if only one type of fermion
        elif not "X" in x_label:
            # we just place underscores after them shits, removing unnecessary ones
            new_label = x_label.split("e")
            for i in range(len(new_label)-1):
                if new_label[i][-1] == "1":
                    new_label[i] = new_label[i][:-1]
            new_label = "e_".join(new_label)
            new_label = new_label.split("h")
            for i in range(len(new_label)-1):
                if new_label[i][-1] == "1":
                    new_label[i] = new_label[i][:-1]
            new_label = "h_".join(new_label)
            return(new_label)
        else:
            #exciton
            # we assume 1 electron, two holes
            cur_occupancies = self.exciton_occupancies[x_label]
            electron_number = cur_occupancies[0][0]
            hole_number = cur_occupancies[1][0]+cur_occupancies[1][1]
            total_charge = hole_number - electron_number
            
            if total_charge == 0:
                charge_label = ""
            elif total_charge == 1:
                charge_label = "+"
            elif total_charge == -1:
                charge_label = "-"
            elif total_charge > 1:
                charge_label = f"{total_charge}+"
            elif total_charge < -1:
                charge_label = f"{-total_charge}-"
            
            if min(electron_number, hole_number) > 1:
                label_prefix = f"{min(electron_number, hole_number)}"
            else:
                label_prefix = ""
            
            new_label = f"{label_prefix}X_{{{cur_occupancies[1][0]},{cur_occupancies[1][1]}}}"
            if charge_label != "":
                new_label += f"^{charge_label}"
            return(new_label)

    
    def output_tikz(self, tikz_string, filename, add_syntax_wrapping = True):
        # prints a tikz string into a file.
        
        os.makedirs("aretdog_outputs", exist_ok=True)
        output_file = open("aretdog_outputs/" + filename + ".tex", "w")
        if add_syntax_wrapping:
            output_file.write("\\documentclass{article}\n\\usepackage{tikz}\n\\begin{document}\n")
        output_file.write(tikz_string)
        if add_syntax_wrapping:
            output_file.write("\n\\end{document}")
        output_file.close()
    
    def tikz_exciton_splitting(self, tikz_picture, exciton_complex, position = (0, 0), scale = 1.0, orientation = 'h',
            line_length = 2, line_separation = 0.5, line_label_distance = 0, line_label_diagonal = 0, line_label_periodicity = 2, exciton_label_distance = 15
            ):
        # This method draws an energy splitting scheme of a excitonic complex (by label) into tikz_picture, an instance of pytikz Picture()
        
        exciton_rep = self.excitons[exciton_complex]
        exciton_rep_reduction = self.separate_constituent_representations(exciton_rep, True)
        
        if orientation in ['h', 'horizontal', 'x']:
            start_x, start_y = position
            cur_x, cur_y = position
            line_index = 0
            for label, rep in exciton_rep_reduction.items():
                proper_label = label[:label.find("[")]
                #print(proper_label)
                if rep.characters["E"] == 1:
                    tikz_picture.draw((cur_x, cur_y), tikz.lineto((cur_x, cur_y+line_length * scale)), tikz.node(f"${tex_readable_irrep_labels[self.name][proper_label]}$", above=f'{line_label_distance-line_label_diagonal*line_index}pt', fill='white', pos=1.0), thick = True)
                elif rep.characters["E"] == 2:
                    tikz_picture.draw((cur_x, cur_y), tikz.lineto((cur_x, cur_y+line_length * scale)), tikz.node(f"${tex_readable_irrep_labels[self.name][proper_label]}$", above=f'{line_label_distance-line_label_diagonal*line_index}pt', fill='white', pos=1.0), ultra_thick = True)
                cur_x += line_separation * scale
                line_index = (line_index + 1) % line_label_periodicity
            cur_x -= line_separation * scale
            tikz_picture.path((start_x, start_y+line_length * scale), tikz.lineto((cur_x, start_y+line_length * scale)), tikz.node(f"${self.tex_readable_exciton_labels(exciton_complex)}$", above=f'{exciton_label_distance}pt', fill='white', pos=0.5, scale=2))
        return((cur_x, start_y+line_length * scale))
    
    def tikz_allowed_transitions(self, tikz_picture, i_X, t_X, i_pos, t_pos, i_max, i_ind, t_max, t_ind, margin=5, spacing = 10, orientation='h'):
        # | [margin] - [box] - [spacing] - [box] - [spacing] - [box] - [margin] |
        # margin and spacing are in percent of box dimension
        
        i_total_length = i_pos[3] - i_pos[1]
        t_total_length = t_pos[3] - t_pos[1]
        
        # b=T/(2*m/100+N+(N-1)*s)
        
        i_box_length = i_total_length / (2.0 * margin / 100.0 + i_max + (i_max - 1) * spacing / 100.0)
        t_box_length = t_total_length / (2.0 * margin / 100.0 + t_max + (t_max - 1) * spacing / 100.0)
        
        i_offset_length = i_box_length * (margin / 100.0 + i_ind * (1.0 + spacing / 100.0))
        t_offset_length = t_box_length * (margin / 100.0 + t_ind * (1.0 + spacing / 100.0))
        
        i_box_pos = [i_pos[0], i_pos[1] + i_offset_length, i_pos[2], i_pos[1] + i_offset_length + i_box_length]
        t_box_pos = [t_pos[0], t_pos[1] + t_offset_length, t_pos[2], t_pos[1] + t_offset_length + t_box_length]
        
        
        # now - because we want all lines parallel, the larger box_length gets limited by the smaller one
        
        if i_box_length > t_box_length:
            i_change = (i_box_length - t_box_length) / 2.0
            i_box_pos[1] += i_change
            i_box_pos[3] -= i_change
            box_length = t_box_length
        elif i_box_length < t_box_length:
            t_change = (t_box_length - i_box_length) / 2.0
            t_box_pos[1] += t_change
            t_box_pos[3] -= t_change
            box_length = i_box_length
        else:
            box_length = i_box_length
        
        #delicious
        i_X_rep = self.excitons[i_X]
        i_X_rep_reduction = self.separate_constituent_representations(i_X_rep, True)
        i_X_irreps = list(i_X_rep_reduction.keys())
        t_X_rep = self.excitons[t_X]
        t_X_rep_reduction = self.separate_constituent_representations(t_X_rep, True)
        t_X_irreps = list(t_X_rep_reduction.keys())
        
        i_E_lvl_N = len(i_X_irreps)
        if i_E_lvl_N > 1:
            i_E_lvl_spacing = (i_pos[2] - i_pos[0]) / (i_E_lvl_N - 1)
        else:
            i_E_lvl_spacing = 0.0
        t_E_lvl_N = len(t_X_irreps)
        if t_E_lvl_N > 1:
            t_E_lvl_spacing = (t_pos[2] - t_pos[0]) / (t_E_lvl_N - 1)
        else:
            t_E_lvl_spacing = 0.0
        
        allowed_transitions = self.allowed_transitions_between_reps(i_X_rep, t_X_rep)
        N_transitions = 0
        for polaris, trans in allowed_transitions.items():
            N_transitions += len(trans[0])
        if N_transitions == 0:
            return(0)
        i_mainline_y = np.linspace(i_box_pos[1], i_box_pos[3], N_transitions)
        t_mainline_y = np.linspace(t_box_pos[1], t_box_pos[3], N_transitions)
        if N_transitions == 1:
            i_mainline_y = np.array([(i_box_pos[3]+i_box_pos[1])/2.0])
            t_mainline_y = np.array([(t_box_pos[3]+t_box_pos[1])/2.0])
        i_mainline_x = (i_box_pos[2]+i_box_pos[0]) / 2.0
        t_mainline_x = (t_box_pos[2]+t_box_pos[0]) / 2.0
        trans_slope = (t_mainline_y[0] - i_mainline_y[0])/(t_mainline_x - i_mainline_x)
        trans_ind = 0
        for polarisation, transitions in allowed_transitions.items():
            if polarisation == "(x,y)":
                linestyle = "solid"
            elif polarisation == "(z)":
                linestyle = "dashed"
            
            for transition in transitions[0]:
                #print(transition, transitions)
                X_labels = transition.split(" -> ")
                i_E_lvl_ind = i_X_irreps.index(X_labels[0])
                t_E_lvl_ind = t_X_irreps.index(X_labels[1])
                
                if i_E_lvl_N > 1:
                    i_trans_pos_x = i_box_pos[0] + i_E_lvl_ind * i_E_lvl_spacing
                else:
                    i_trans_pos_x = i_mainline_x
                i_trans_pos_y = i_mainline_y[trans_ind] + (i_trans_pos_x - i_mainline_x) * trans_slope
                if t_E_lvl_N > 1:
                    t_trans_pos_x = t_box_pos[0] + t_E_lvl_ind * t_E_lvl_spacing
                else:
                    t_trans_pos_x = t_mainline_x
                t_trans_pos_y = t_mainline_y[trans_ind] + (t_trans_pos_x - t_mainline_x) * trans_slope
                
                tikz_picture.draw((i_trans_pos_x, i_trans_pos_y), tikz.lineto((t_trans_pos_x, t_trans_pos_y)), tikz.node(f"$\\nu_{trans_ind+1}$", fill='white', pos=0.5), thick = True, opt=linestyle + ",->")
                
                trans_ind += 1
            
    
    def tikz_decay_diagram_print(self, exciton_complex,
            line_length = 2, line_separation = 0.5, line_label_distance = 0, line_label_diagonal = 0, line_label_periodicity = 2, exciton_label_distance = 15):
        
        decay_diagram_pic = tikz.Picture(scale=1)
        
        environment_width = 8
        environment_height = 6
        environment_pos = (0, 0)
        environment_pos_x, environment_pos_y = environment_pos
        base_line_length = 2
        #complex_distance = 6 #distance from A to B in A->B
        #complex_separation = 4 #perpendicular distance between alternative complexes
        
        encompassed_complexes = [[exciton_complex]] #[index in diagram:[complex 1, complex 2, ...]]
        
        number_of_transitions_in = [[0]] # same shape
        number_of_transitions_in_taken = [[0]]
        
        max_complexes_in_row = 1
        
        is_this_the_end = False
        while(not is_this_the_end):
            is_this_the_end = True
            new_encompassed_complexes = []
            new_number_of_transitions_in = []
            new_number_of_transitions_in_taken = []
            for X in encompassed_complexes[-1]:
                for new_X in self.transition_chain[X]:
                    if not new_X in new_encompassed_complexes:
                        new_encompassed_complexes.append(new_X)
                        new_number_of_transitions_in.append(1)
                        new_number_of_transitions_in_taken.append(0)
                    else:
                        new_number_of_transitions_in[new_encompassed_complexes.index(new_X)] += 1
            if len(new_encompassed_complexes) > 0:
                is_this_the_end = False
                encompassed_complexes.append(new_encompassed_complexes)
                number_of_transitions_in.append(new_number_of_transitions_in)
                number_of_transitions_in_taken.append(new_number_of_transitions_in_taken)
                if len(new_encompassed_complexes) > max_complexes_in_row:
                    max_complexes_in_row = len(new_encompassed_complexes)
        
        complex_distance = environment_width / (len(number_of_transitions_in) - 1)
        
        position_of_complexes = []
        
        for i in range(len(encompassed_complexes)):
            position_of_complexes.append([])
            N = len(encompassed_complexes[i])
            cur_pos_x = environment_pos_x + i * complex_distance
            
            pos_y_list = np.linspace(environment_pos_y, environment_pos_y + environment_height, len(encompassed_complexes[i]))
            if len(encompassed_complexes[i]) == 1:
                pos_y_list = np.array([environment_pos_y / 2.0 + environment_height / 2.0])
            
            cur_line_length = base_line_length * max_complexes_in_row / len(encompassed_complexes[i])
            for j in range(N):
                cur_pos_y = pos_y_list[j]-cur_line_length/2.0#(j-N/2.0) * complex_separation
                final_pos_x, final_pos_y = self.tikz_exciton_splitting(decay_diagram_pic, encompassed_complexes[i][j], (cur_pos_x, cur_pos_y),scale=1.0, orientation='h', line_length = cur_line_length)
                position_of_complexes[-1].append([cur_pos_x, cur_pos_y, final_pos_x, final_pos_y])
        
        # now we add the transition lines
        for i in range(len(encompassed_complexes)-1):
            for j in range(len(encompassed_complexes[i])):
                cur_exciton_complex = encompassed_complexes[i][j]
                number_of_transitions_out = len(self.transition_chain[cur_exciton_complex])
                for k in range(number_of_transitions_out):
                    target_exciton = self.transition_chain[cur_exciton_complex][k]
                    target_exciton_list_index = encompassed_complexes[i+1].index(target_exciton)
                    self.tikz_allowed_transitions(decay_diagram_pic, cur_exciton_complex, target_exciton, position_of_complexes[i][j], position_of_complexes[i+1][target_exciton_list_index], number_of_transitions_out, k, number_of_transitions_in[i+1][target_exciton_list_index], number_of_transitions_in_taken[i+1][target_exciton_list_index], margin=15, spacing = 20)
                    number_of_transitions_in_taken[i+1][target_exciton_list_index] += 1
                #for target_complex in self.transition_chain[encompassed_complexes[i][j]]
        
        self.output_tikz(decay_diagram_pic.code(), "decay_diagram_in_" + self.name)
                
            
            
        
    
        
        
        
    
    
    
