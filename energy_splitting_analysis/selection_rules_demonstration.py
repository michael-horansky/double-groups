


from aretdog.class_group import *




QDG = Group("D3h QD")

QDG.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "m" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True), "m'" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], True)})

QDG.print_character_table()

#print("gamma(j=1/2) x gamma(j=1/2) rep decomp:", QDG.reduce_representation(QDG.angular_representation(1/2)*QDG.angular_representation(1/2))[1])
#print("gamma(j=1/2) x gamma(j=5/2) rep decomp:", QDG.reduce_representation(QDG.angular_representation(1/2)*QDG.irrep_characters["E_3(j=5/2)"])[1])



#two_X_11 = QDG.irrep_characters["E_1(j=1/2)"] * (QDG.irrep_characters["1E(j=3/2)"] + QDG.irrep_characters["2E(j=3/2)"])
#X_10 = QDG.irrep_characters["E_1(j=1/2)"] * (QDG.irrep_characters["1E(j=3/2)"] + QDG.irrep_characters["2E(j=3/2)"])
X_01 = QDG.irrep_characters["E_2(j=1/2)"] * QDG.irrep_characters["E_3(j=5/2)"]
vacuum = QDG.irrep_characters["A_1"]

print("X_01:", QDG.reduce_representation(X_01)[1])
print("vacuum:", QDG.reduce_representation(vacuum)[1])


transitions = QDG.allowed_transitions_between_reps(X_01, vacuum)

for polarisation in transitions:
    print(polarisation + "-polarised transitions in X_01 -> vacuum:")
    print("  Allowed: ", transitions[polarisation][0])
    print("  Dark:    ", transitions[polarisation][1])
