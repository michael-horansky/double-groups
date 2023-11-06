


from aretdog.class_group import *

C_i = Group("C_i")
C_i.generate_group({"E" : E, "i" : ImproperRotation([1.0, 0.0, 0.0], [0, 1], True)})
C_i.print_character_table()


D2_group = Group("D2")

D2_group.generate_double_group({"E" : E, "a" : Cz_2, "b" : Cy_2})


D2_group.print_character_table()


T_group = Group("T")
T_group.generate_double_group({"E" : E, "Cz_2" : Cz_2, "Cy_2" : Cy_2, "C'_3" : ImproperRotation([1.0, 1.0, 1.0], [1, 3], False)})


T_group.print_character_table()

print("j = 1/2:",T_group.angular_representation(1/2), T_group.reduce_representation(T_group.angular_representation(1/2)))
print("j = 3/2:", T_group.reduce_representation(T_group.angular_representation(3/2)))
print("j = 5/2:", T_group.reduce_representation(T_group.angular_representation(5/2)))
print("j = 7/2:", T_group.reduce_representation(T_group.angular_representation(7/2)))

T_group.add_subgroup(D2_group)

print("---------------------- symmetry breakage -------------------")

og_rep = T_group.irrep_characters["T_1"]
rep = T_group.rep_to_subgroup_rep("D2", og_rep)

print(og_rep)
print(rep)
print(D2_group.reduce_representation(rep))


"""
Th_group = Group("Th")
Th_group.generate_double_group({"E" : E, "Cz_2" : Cz_2, "Cy_2" : Cy_2, "C'_3" : ImproperRotation([1.0, 1.0, 1.0], [1, 3], False), "i" : ImproperRotation([1.0, 0.0, 0.0], [0.0, 1.0], True)})
Th_group.print_character_table()"""



print("------------------ C3v analysis ------------------------")

C3v_group = Group("C3v")
C3v_group.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C3v_group.print_character_table()




D3h_group = Group("D3h")
D3h_group.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "m" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True), "m'" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], True)})
D3h_group.print_character_table()


D3h_group.add_subgroup(C3v_group)

two_X_11 = C3v_group.irrep_characters["E_1(j=1/2)"] * (C3v_group.irrep_characters["1E(j=3/2)"] + C3v_group.irrep_characters["2E(j=3/2)"])
X_10 = C3v_group.irrep_characters["E_1(j=1/2)"] * (C3v_group.irrep_characters["1E(j=3/2)"] + C3v_group.irrep_characters["2E(j=3/2)"])
X_01 = C3v_group.irrep_characters["E_1(j=1/2)"] * C3v_group.irrep_characters["E_1(j=1/2)"]
vacuum = C3v_group.irrep_characters["A_1"]
C3v_dipole_rep = C3v_group.irrep_characters["E_1"] # Karlsson p 13 - in an optic transition, we do dipole approximation and Wigner-Eckart (at least for isotropic polarisation)

print("2X_11 has energy leveles as", C3v_group.reduce_representation(two_X_11)[1])
print("X_10 has energy leveles as", C3v_group.reduce_representation(X_10)[1])
print("X_01 has energy leveles as", C3v_group.reduce_representation(X_01)[1])

X_11_to_X_10 = C3v_group.allowed_transitions_between_reps(two_X_11, X_10)

for polarisation in X_11_to_X_10:
    print(polarisation + "-polarised transitions in 2X_11 -> X_10:")
    print("  Allowed: ", X_11_to_X_10[polarisation][0])
    print("  Dark:    ", X_11_to_X_10[polarisation][1])

X_01_to_vacuum = C3v_group.allowed_transitions_between_reps(X_01, vacuum)
for polarisation in X_01_to_vacuum:
    print(polarisation + "-polarised transitions in X_01 -> vacuum:")
    print("  Allowed: ", X_01_to_vacuum[polarisation][0])
    print("  Dark:    ", X_01_to_vacuum[polarisation][1])

"""
print("Allowed isotropically-polarised transitions in 2X_11 -> X_10:", allowed)
print("Dark isotropically-polarised transitions in 2X_11 -> X_10:", dark)
allowed, dark = C3v_group.allowed_transitions_between_reps(X_01, vacuum, C3v_dipole_rep)
print("Allowed isotropically-polarised transitions in 2X_11 -> X_10:", allowed)
print("Dark isotropically-polarised transitions in 2X_11 -> X_10:", dark)"""



