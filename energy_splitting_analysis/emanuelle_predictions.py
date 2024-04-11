


from aretdog.class_QD_group import *



print("Constructing the double group C3v...")
C3v = QDGroup("C3v QD")
C3v.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C3v.print_character_table()
C3v.classify_excitons(max_e = 3, max_h = 3)
#C3v.print_exciton_complexes()
C3v.find_transition_chain()
#print(C3v.transition_chain)

print("Constructing the double group C6v...")
C6v = QDGroup("C6v QD")
C6v.generate_double_group({"E" : E, "Cz_6" : ImproperRotation([0.0, 0.0, 1.0], [1, 6], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C6v.print_character_table()
C6v.classify_excitons(max_e = 3, max_h = 3)
#C6v.print_exciton_complexes()
C6v.find_transition_chain()
#print(C3v.transition_chain)
C6v.add_subgroup(C3v)

print("Constructing the double group D3h...")
D3h = QDGroup("D3h QD")
D3h.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "m" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True), "m'" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], True)})
D3h.print_character_table()
D3h.classify_excitons(max_e = 3, max_h = 3)
#D3h.print_exciton_complexes()

print("Constructing the double group D3d...")
D3d = QDGroup("D3d QD")
D3d.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "C_2" : ImproperRotation([1.0, -1.0, 0.0], [1, 2], False), "i" : ImproperRotation([0.0, 0.0, 1.0], [1, 1], True)})
D3d.print_character_table()
D3d.classify_excitons(max_e = 3, max_h = 3)
#D3d.print_exciton_complexes()

print("Constructing the double group Oh...")
Oh = QDGroup("Oh QD")
Oh.generate_double_group({"E" : E,
                          "Cz_2" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], False),
                          "Cy_2" : ImproperRotation([0.0, 1.0, 0.0], [1, 2], False),
                          "C'_2" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], False),
                          "C_3"   : ImproperRotation([1.0, 1.0, 1.0], [1, 3], False),
                          "sigma" : ImproperRotation([1.0, 0.0, 0.0], [1, 1], True)})
Oh.print_character_table()
Oh.classify_excitons(max_e = 3, max_h = 3)
Oh.print_exciton_complexes()

print("Constructing the double group Td...")
Td = QDGroup("Td QD")
Td.generate_double_group({"E" : E,
                          "Cz_2" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], False),
                          "Cy_2" : ImproperRotation([0.0, 1.0, 0.0], [1, 2], False),
                          "C_3"   : ImproperRotation([1.0, 1.0, 1.0], [1, 3], False),
                          "m" : ImproperRotation([1.0, -1.0, 0.0], [1, 2], True)})
Td.print_character_table()
Td.classify_excitons(max_e = 3, max_h = 3)
Td.print_exciton_complexes()



# Something doesnt add up
# h1 + h2 = j=3/2 wigner rep
"""
print(" ------------ wigner check --------------")
wigner_rep = C3v.find_wigner_representation(3/2)
print(C3v.reduce_representation(C3v.angular_representation(3/2) * C3v.irrep_characters[C3v.inversion_irrep])[1])
wigner_rep = C6v.find_wigner_representation(3/2)

print(C6v.reduce_representation(C6v.irrep_characters["E_1(j=1/2)"] * C6v.irrep_characters[C6v.inversion_irrep])[1])
print(C6v.reduce_representation(C6v.irrep_characters["E_3(j=3/2)"] * C6v.irrep_characters[C6v.inversion_irrep])[1])

print(C6v.reduce_representation(C6v.angular_representation(3/2) * C6v.irrep_characters[C6v.inversion_irrep])[1])
print(D3h.reduce_representation(D3h.angular_representation(3/2) * D3h.irrep_characters[D3h.inversion_irrep])[1])"""

print(" ------------ differing predictions ---------")
# both transition chains are equal - we will now compare element by element, and see which ones differ in the number of emission lines
inequal_transitions = [] # each element is a two-element list [init state, final state]
for initial_state, final_states in C3v.transition_chain.items():
    for final_state in final_states:
        #print(initial_state, final_state)
        C3v_transitions = C3v.allowed_transitions_between_reps(C3v.excitons[initial_state], C3v.excitons[final_state])
        C6v_transitions = C6v.allowed_transitions_between_reps(C6v.excitons[initial_state], C6v.excitons[final_state])
        #print(C3v_transitions)
        equal_transitions = True
        for polarisation in C3v_transitions.keys():
            #print(len(C3v_transitions[polarisation][0]), len(C6v_transitions[polarisation][0]))
            if len(C3v_transitions[polarisation][0]) != len(C6v_transitions[polarisation][0]):
                equal_transitions = False
                break
        if not equal_transitions:
            inequal_transitions.append([initial_state, final_state])

print("Inequal transitions for C3v -> C6v symmetry elevation")
for inequal_transition in inequal_transitions:
    i_state = inequal_transition[0]
    f_state = inequal_transition[1]
    print(f"  {i_state} -> {f_state}:")
    C3v_transitions = C3v.allowed_transitions_between_reps(C3v.excitons[i_state], C3v.excitons[f_state])
    C6v_transitions = C6v.allowed_transitions_between_reps(C6v.excitons[i_state], C6v.excitons[f_state])
    for polarisation in C3v_transitions.keys():
        if len(C3v_transitions[polarisation][0]) != len(C6v_transitions[polarisation][0]):
            print(f"    {polarisation}-polar.: C3v = {len(C3v_transitions[polarisation][0])}; C6v = {len(C6v_transitions[polarisation][0])}")
        #print(f"      C3v: {'; '.join(C3v_transitions[polarisation][0])}")
        #print(f"      C6v: {'; '.join(C6v_transitions[polarisation][0])}")

#print_table("test", [1, 56, 64321], ["a", 72], [[1, 2, 3], [False, 5, 6]], [1])

def polarisation_overlap(beam_polarisation, detector_polarisation):
    # both in the notation "a,b,c"
    # two different circular polarisations have zero overlap because the beams irradiate in different directions
    if beam_polarisation == "x,y,z":
        if detector_polarisation in ["x", "y", "z"]:
            return(0.5)
        else:
            return(1)
    elif beam_polarisation in ["x,y", "x,z", "y,z"]:
        if detector_polarisation in ["x", "y", "z"]:
            if detector_polarisation in beam_polarisation:
                return(0.5)
            else:
                return(0)
        else:
            if detector_polarisation == beam_polarisation:
                return(1)
            else:
                return(0) #because of detector orientation - we assume the beams to be localised. Check this.
    elif beam_polarisation in ["x", "y", "z"]:
        if detector_polarisation in ["x", "y", "z"]:
            if detector_polarisation == beam_polarisation:
                return(1)
            else:
                return(0)
        else:
            if beam_polarisation in detector_polarisation:
                return(1)
            else:
                return(0)

def OLD_compare_emission_line_numbers(list_of_groups, detector_polarisations = ["x,y", "z"]):
    # all groups must be initialised to the same decay chain
    # all groups must polarise light in the same form of division of basis partners
    # TODO not anymore - instead we calculate the overlap with the set of detectors, with the coefficient included in the table data
    # two unequal linearly polarised beams have 0 overlap
    # unpolarised light and a linearly polarised beam has 0.5 overlap (Malus's law)
    
    # here we obtain the shared polarisations
    raw_polarisations = list(list_of_groups[0].cartesian_basis.values())
    polarisations = []
    for same_irrep_polaris in raw_polarisations:
        for polaris in same_irrep_polaris:
            polarisations.append(polaris)
    for i in range(len(polarisations)):
        polarisations[i] = "(" + polarisations[i] + ")"
    
    headers = []
    subtable_borders = []
    for i in range(len(polarisations)):
        for qd_group in list_of_groups:
            headers.append(f"{qd_group.name} {polarisations[i]}")
        subtable_borders.append((i+1) * len(list_of_groups) - 1)
    subtable_borders = subtable_borders[:-1]
    
    row_labels = []
    table_data = []
    for i_state, f_states in list_of_groups[0].transition_chain.items():
        for f_state in f_states:
            row_labels.append(f"{i_state} -> {f_state}:")
            table_data.append([])
            cur_group_transitions = []
            for qd_group in list_of_groups:
                cur_group_transitions.append(qd_group.allowed_transitions_between_reps(qd_group.excitons[i_state], qd_group.excitons[f_state]))
            for polarisation in polarisations:
                for i_group in range(len(list_of_groups)):
                    table_data[-1].append(len(cur_group_transitions[i_group][polarisation][0]))
    
    print_table("# of em. lines", headers, row_labels, table_data, subtable_borders)

def compare_emission_line_numbers(list_of_groups, detector_polarisations = ["x,y", "z"]):
    # all groups must be initialised to the same decay chain not true anymore
    # all groups must polarise light in the same form of division of basis partners not true anymore lessgoo
    # TODO not anymore - instead we calculate the overlap with the set of detectors, with the coefficient included in the table data
    # two unequal linearly polarised beams have 0 overlap
    # unpolarised light and a linearly polarised beam has 0.5 overlap (Malus's law)
    
    # here we obtain the shared polarisations
    """raw_polarisations = list(list_of_groups[0].cartesian_basis.values())
    polarisations = []
    for same_irrep_polaris in raw_polarisations:
        for polaris in same_irrep_polaris:
            polarisations.append(polaris)
    for i in range(len(polarisations)):
        polarisations[i] = "(" + polarisations[i] + ") """
    
    headers = []
    subtable_borders = []
    for i in range(len(detector_polarisations)):
        for qd_group in list_of_groups:
            headers.append(f"{qd_group.name} ({detector_polarisations[i]})")
        subtable_borders.append((i+1) * len(list_of_groups) - 1)
    subtable_borders = subtable_borders[:-1]
    
    row_labels = []
    table_data = []
    for i_state, f_states in list_of_groups[0].transition_chain.items():
        for f_state in f_states:
            row_labels.append(f"{i_state} -> {f_state}:")
            table_data.append([])
            cur_group_transitions = []
            for qd_group in list_of_groups:
                # Here we calculate the reduction to states for groups with fewer singular fermionic characters
                # first, if numbers of characters agree, we just use the label as the dict key
                e_char_ref = len(list_of_groups[0].electrons)
                h_char_ref = len(list_of_groups[0].holes)
                e_char_cur = len(qd_group.electrons)
                h_char_cur = len(qd_group.holes)
                
                if e_char_ref == e_char_cur and h_char_ref == h_char_cur:
                    cur_group_transitions.append(qd_group.allowed_transitions_between_reps(qd_group.excitons[i_state], qd_group.excitons[f_state]))
                else:
                    if e_char_cur == 1 and h_char_cur == 1:
                        i_e_occupancy = sum(list_of_groups[0].exciton_occupancies[i_state][0])
                        i_h_occupancy = sum(list_of_groups[0].exciton_occupancies[i_state][1])
                        f_e_occupancy = sum(list_of_groups[0].exciton_occupancies[f_state][0])
                        f_h_occupancy = sum(list_of_groups[0].exciton_occupancies[f_state][1])
                        new_i_state = exciton_label_from_occupancies([i_e_occupancy], [i_h_occupancy])
                        new_f_state = exciton_label_from_occupancies([f_e_occupancy], [f_h_occupancy])
                        cur_group_transitions.append(qd_group.allowed_transitions_between_reps(qd_group.excitons[new_i_state], qd_group.excitons[new_f_state]))
            
            for d_polarisation in detector_polarisations:
                for i_group in range(len(list_of_groups)):
                    final_cell_data = ""
                    number_of_lines_full_intensity = 0
                    number_of_lines_half_intensity = 0
                    for b_polarisation in cur_group_transitions[i_group].keys():
                        number_of_lines = len(cur_group_transitions[i_group][b_polarisation][0])
                        intensity_of_lines = polarisation_overlap(b_polarisation[1:-1], d_polarisation)
                        if intensity_of_lines == 0 or number_of_lines == 0:
                            continue
                        elif intensity_of_lines == 1:
                            number_of_lines_full_intensity += number_of_lines
                        elif intensity_of_lines == 0.5:
                            number_of_lines_half_intensity += number_of_lines
                        else:
                            print("Something strange happened.")
                    
                    if number_of_lines_half_intensity == 0:
                        if number_of_lines_full_intensity == 0:
                            final_cell_data = "0"
                        else:
                            final_cell_data = str(number_of_lines_full_intensity)
                    else:
                        if number_of_lines_full_intensity == 0:
                            final_cell_data = str(number_of_lines_half_intensity) + " (i.1/2)"
                        else:
                            final_cell_data = f"{number_of_lines_full_intensity} + {number_of_lines_half_intensity} (i.1/2)"
                    table_data[-1].append(final_cell_data)
    
    print_table("# of em. lines", headers, row_labels, table_data, subtable_borders)

print("---------------------------------------------")
#OLD_compare_emission_line_numbers([C3v, C6v, D3h, D3d])
#compare_emission_line_numbers([C3v, C6v, D3h, D3d])
compare_emission_line_numbers([C3v, C6v, Oh, Td])

#for i in range(len(Oh.multihole_states)):
#    print(f"{i+1} holes in Oh would transform as {Oh.reduce_representation(Oh.multihole_states[i])[1]}")



#lel = Oh.angular_representation(2, "g") + Oh.angular_representation(0, "g")
#print(Oh.reduce_representation(lel)[1])

#D3h.output_tikz("\\draw[gray, thick] (-1,2) -- (2,-4);\n\\draw[gray, thick] (-1,-1) -- (2,2);\n\\filldraw[black] (0,0) circle (2pt) node[anchor=west]{Intersection point};", "test", True)
#C3v.tikz_decay_diagram_print("2X[2][2,1]")
#C6v.tikz_decay_diagram_print("2X[2][2,1]")
#D3h.tikz_decay_diagram_print("2X[3][1,1]")

"""
print("---------------------- debug ---------------------")

d3h_e5 = D3h.irrep_characters["E_3(j=5/2)"]
d3h_e1 = D3h.irrep_characters["E_2(j=1/2)"]
d3h_xy = D3h.irrep_characters["E_2"]

trans_rep = d3h_e1 * d3h_e5 * d3h_xy
print(D3h.reduce_representation(trans_rep), D3h.does_rep_contain_identity(trans_rep))
id_coef = 0
        
for i in range(len(D3h.conjugacy_class_names)):
    id_coef += trans_rep.characters[D3h.conjugacy_class_names[i]] * D3h.conjugacy_class_sizes[i]
id_coef /= D3h.order

print("id coef =", id_coef)

lol = D3h.allowed_transitions_between_reps(d3h_e5, d3h_e1)
for polaris, trans in lol.items():
    print(f"  {polaris}-polar.: {'; '.join(trans[0])}")"""
