


from aretdog.class_QD_group import *

#We orientate the C3v rotation axis along [111]

C3v = QDGroup("C3v QD")
C3v.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([1.0, 1.0, 1.0], [1, 3], False), "sigma" : ImproperRotation([-1.0, 1.0, 0.0], [1, 2], True)})



Oh = QDGroup("Oh QD")
Oh.generate_double_group({"E" : E,
                          "Cz_2" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], False),
                          "Cy_2" : ImproperRotation([0.0, 1.0, 0.0], [1, 2], False),
                          "C'_2" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], False),
                          "C_3"   : ImproperRotation([1.0, 1.0, 1.0], [1, 3], False),
                          "sigma" : ImproperRotation([1.0, 0.0, 0.0], [1, 1], True)})
Oh.print_character_table()

Td = QDGroup("Td QD")
Td.generate_double_group({"E" : E,
                          "Cz_2" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], False),
                          "Cy_2" : ImproperRotation([0.0, 1.0, 0.0], [1, 2], False),
                          "C_3"   : ImproperRotation([1.0, 1.0, 1.0], [1, 3], False),
                          "m" : ImproperRotation([1.0, -1.0, 0.0], [1, 2], True)})
Td.print_character_table()

# m_110 = 2_110 * i, hence in this (and only this) orientation, every symmetry transformation in C_3v is also present in O_h
