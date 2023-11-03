


from aretdog.groups import *



D2_group = Group("D2")

D2_group.generate_double_group({"E" : E, "Cz_2" : Cz_2, "Cy_2" : Cy_2})


D2_group.print_character_table()

T_group = Group("T")
T_group.generate_double_group({"E" : E, "Cz_2" : Cz_2, "Cy_2" : Cy_2, "C'_3" : ImproperRotation([1.0, 1.0, 1.0], [1, 3], False)})
#print(T_group.conjugacy_class_names)
#print(T_group.irrep_names)
#print(T_group.character_table)

T_group.print_character_table()


a = T_group.angular_representation(5/2)
print(T_group.reduce_representation(a))

T_group.add_subgroup(D2_group)

print(T_group.subgroup_conjugacy_relations)

"""
print("---------------------- symmetry breakage -------------------")

og_rep = T_group.character_table[-1]
rep = T_group.rep_to_subgroup_rep("D2", og_rep)
print("In T:", og_rep)
print("In D2:", rep)
print(D2_group.reduce_representation(rep))


print("------------------ improper groups study -------------------")

Ci_group = Group("C_i")

Ci_group.generate_double_group({"E" : E, "i" : ImproperRotation([0.0, 0.0, 1.0], [0, 1], True)})

#Ci_group.initialize_from_multiplication_table()

Ci_group.print_character_table()

Ci_group = Group("C_s")
Ci_group.generate_double_group({"E" : E, "i" : ImproperRotation([1.0, 0.0, 0.0], [1, 2], True)})
Ci_group.print_character_table()

C3v_group = Group("C3v")
C3v_group.generate_double_group({"E" : E, "C_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C3v_group.print_character_table()"""
