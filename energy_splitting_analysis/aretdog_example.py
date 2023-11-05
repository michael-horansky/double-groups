


from aretdog.class_group import *




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

allowed, dark = C3v_group.allowed_transitions_between_reps(two_X_11, X_10, C3v_dipole_rep)
print("Allowed isotropically-polarised transitions in 2X_11 -> X_10:", allowed)
print("Dark isotropically-polarised transitions in 2X_11 -> X_10:", dark)
allowed, dark = C3v_group.allowed_transitions_between_reps(X_01, vacuum, C3v_dipole_rep)
print("Allowed isotropically-polarised transitions in 2X_11 -> X_10:", allowed)
print("Dark isotropically-polarised transitions in 2X_11 -> X_10:", dark)

E_12 = D3h_group.character_table[7] #E1/2
E_52 = D3h_group.character_table[8] #E5/2

E_12_in_C3v = D3h_group.rep_to_subgroup_rep("C3v", E_12)
E_52_in_C3v = D3h_group.rep_to_subgroup_rep("C3v", E_52)

