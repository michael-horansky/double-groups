


from aretdog.class_QD_group import *




C3v = QDGroup("C3v QD")
C3v.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C3v.print_character_table()
C3v.classify_excitons(max_e = 3, max_h = 3)
C3v.print_exciton_complexes()
C3v.find_transition_chain()
print(C3v.transition_chain)

C6v = QDGroup("C6v QD")
C6v.generate_double_group({"E" : E, "Cz_6" : ImproperRotation([0.0, 0.0, 1.0], [1, 6], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C6v.print_character_table()
C6v.classify_excitons(max_e = 3, max_h = 3)
C6v.print_exciton_complexes()
C6v.find_transition_chain()
print(C3v.transition_chain)
C6v.add_subgroup(C3v)

D3h = QDGroup("D3h QD")
D3h.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "m" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True), "m'" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], True)})
#D3h.print_character_table()
D3h.classify_excitons(max_e = 2, max_h = 2)
#D3h.print_exciton_complexes()

# Something doesnt add up
# h1 + h2 = j=3/2 wigner rep
print(" ------------ wigner check --------------")
wigner_rep = C3v.find_wigner_representation(3/2)
print(C3v.reduce_representation(C3v.angular_representation(3/2) * C3v.irrep_characters[C3v.inversion_irrep])[1])
wigner_rep = C6v.find_wigner_representation(3/2)

print(C6v.reduce_representation(C6v.irrep_characters["E_1(j=1/2)"] * C6v.irrep_characters[C6v.inversion_irrep])[1])
print(C6v.reduce_representation(C6v.irrep_characters["E_3(j=3/2)"] * C6v.irrep_characters[C6v.inversion_irrep])[1])

print(C6v.reduce_representation(C6v.angular_representation(3/2) * C6v.irrep_characters[C6v.inversion_irrep])[1])
print(D3h.reduce_representation(D3h.angular_representation(3/2) * D3h.irrep_characters[D3h.inversion_irrep])[1])

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
