


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



Th_group = Group("Th")
Th_group.generate_double_group({"E" : E, "Cz_2" : Cz_2, "Cy_2" : Cy_2, "C'_3" : ImproperRotation([1.0, 1.0, 1.0], [1, 3], False), "i" : ImproperRotation([1.0, 0.0, 0.0], [0.0, 1.0], True)})
Th_group.print_character_table()



print("------------------ C3v analysis ------------------------")

C3v_group = Group("C3v")
C3v_group.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C3v_group.print_character_table()




D3h_group = Group("D3h")
D3h_group.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "m" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True), "m'" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], True)})
D3h_group.print_character_table()

D3h_group.add_subgroup(C3v_group)

E_12 = D3h_group.character_table[7] #E1/2
E_52 = D3h_group.character_table[8] #E5/2

E_12_in_C3v = D3h_group.rep_to_subgroup_rep("C3v", E_12)
E_52_in_C3v = D3h_group.rep_to_subgroup_rep("C3v", E_52)

print(C3v_group.reduce_representation(E_12_in_C3v))
print(C3v_group.reduce_representation(E_52_in_C3v)) #this is correct!!


#for g in D3h_group.group_elements:
#    print(D3h_group.element_spatial_properties[g][1])



"""
a = T_group.angular_representation(5/2)
print(T_group.reduce_representation(a))


print(T_group.subgroup_conjugacy_relations)
"""
"""



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
