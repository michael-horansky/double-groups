


from aretdog.class_QD_group import *




C3v = QDGroup("C3v QD")
C3v.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C3v.print_character_table()
C3v.classify_excitons(max_e = 3, max_h = 2)
C3v.print_exciton_complexes()

D3h = QDGroup("D3h QD")
D3h.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "m" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True), "m'" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], True)})
D3h.print_character_table()
D3h.classify_excitons(max_e = 3, max_h = 2)
D3h.print_exciton_complexes()

D3h.find_transition_chain()
print(D3h.transition_chain)

"""

print("X_01:", QDG.reduce_representation(X_01)[1])
print("vacuum:", QDG.reduce_representation(vacuum)[1])


transitions = QDG.allowed_transitions_between_reps(X_01, vacuum)

for polarisation in transitions:
    print(polarisation + "-polarised transitions in X_01 -> vacuum:")
    print("  Allowed: ", transitions[polarisation][0])
    print("  Dark:    ", transitions[polarisation][1])
"""
