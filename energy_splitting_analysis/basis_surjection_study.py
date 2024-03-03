


from aretdog.class_QD_group import *




C3v = QDGroup("C3v QD")
C3v.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C3v.print_character_table()
C3v.classify_excitons(max_e = 2, max_h = 3)
C3v.print_exciton_complexes()
#C3v.find_transition_chain()
#print(C3v.transition_chain)

C6v = QDGroup("C6v QD")
C6v.generate_double_group({"E" : E, "Cz_6" : ImproperRotation([0.0, 0.0, 1.0], [1, 6], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C6v.print_character_table()
C6v.classify_excitons(max_e = 2, max_h = 3)
C6v.print_exciton_complexes()
#C3v.find_transition_chain()
#print(C3v.transition_chain)
C6v.add_subgroup(C3v)

D3h = QDGroup("D3h QD")
D3h.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "m" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True), "m'" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], True)})
D3h.print_character_table()
D3h.classify_excitons(max_e = 2, max_h = 3)
D3h.print_exciton_complexes()
#D3h.find_transition_chain()
#print(D3h.transition_chain)
D3h.add_subgroup(C3v)


print("----------------- Basis surjection study -------------------")

lmao = C6v.basis_surjection_to_subgroup("C3v QD", True)
print(lmao)

print("From this we see that for X_1,1 one of the E_1/2 elevates to E_5/2, which gets suppressed in the transition from 2X_2,1 (which is E_1/2 elevating to E_5/2)")

#C6v.tikz_basis_surjection_diagram("C3v QD")
D3h.tikz_basis_surjection_diagram("C3v QD", position = (15, 0))

