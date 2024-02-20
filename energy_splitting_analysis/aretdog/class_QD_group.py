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


from tikz import *
import os
from .class_group import *


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
    
    def output_tikz(self, tikz_string, filename, add_syntax_wrapping = True):
        # prints a tikz string into a file.
        
        os.makedirs("aretdog_outputs", exist_ok=True)
        output_file = open("aretdog_outputs/" + filename + ".tex", "w")
        if add_syntax_wrapping:
            output_file.write("\\documentclass{article}\n\\usepackage{tikz}\n\\begin{document}\n\\begin{tikzpicture}\n")
        output_file.write(tikz_string)
        if add_syntax_wrapping:
            output_file.write("\n\\end{tikzpicture}\n\\end{document}")
        output_file.close()
        
    
        
        
        
    
    
    
