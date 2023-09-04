from improper_rotations import *


#TODO add basis functions idk


class Group():
    
    def __init__(self, conjugacy_classes, representations, character_table):
        # conjugacy_classes and representations are lists of labels for these respectively, OR integer lengths of these if no labels are necessary
        self.conjugacy_classes = conjugacy_classes
        if type(conjugacy_classes) == int:
            self.conjugacy_classes = []
            for i in range(conjugacy_classes):
                self.conjugacy_classes.append(str(i))
        
        self.representations = representations
        if type(representations) == int:
            self.representations = []
            for i in range(representations):
                self.representations.append(str(i))
        
        # Character table MUST be the dimension [len(representations)][len(conjugacy_classes)]
        
        self.character_table = np.array(character_table)
        
        self.conjugacy_class_elements = [] #i-th element is the list of symmetry operations in the i-th CC
        for i in range(len(self.conjugacy_classes)):
            self.conjugacy_class_elements.append([])
        
        self.group_elements = {} # a dictionary {symmetry operation : conjugacy class index}
        
    
    def add_symmetry_operation(self, group_element, conjugacy_class_index):
        self.conjugacy_class_elements[conjugacy_class_index].append(group_element)
        self.group_elements[group_element] = conjugacy_class_index
        
        


# ------------------ Example: D_4 -> D_2 -----------------

E = ImproperRotation([1.0, 0.0, 0.0], 0.0, False)
Cz_4 = ImproperRotation([0.0, 0.0, 1.0], np.pi / 2.0, False)
Cz_2 = ImproperRotation([0.0, 0.0, 1.0], np.pi, False)
Cz_43 = ImproperRotation([0.0, 0.0, 1.0], 3.0 * np.pi / 2.0, False)

Cx_2 = ImproperRotation([1.0, 0.0, 0.0], np.pi, False)
Cy_2 = ImproperRotation([0.0, 1.0, 0.0], np.pi, False)
Cdia_1 = Cz_4 + Cy_2
Cdia_2 = Cy_2 + Cz_4

print(Cdia_2 == Cz_4 + Cz_4 + Cz_4 + Cy_2) #LETSGOOO THIS IS TRUEEEEEEE




"""
x = ImproperRotation([0, 0, 1], np.pi, False)
y = ImproperRotation([0, 1, 0], np.pi, True)

print("x =", x.cartesian_rep())
print("y =", y.cartesian_rep())
z = x + y
print("x + y =", z.cartesian_rep())
print("x + y axis =", z.axis, "; x + y angle =", z.angle)

w = ImproperRotation([1, 0, 0], np.pi, False)

print(x == w)
print(z == w)

print("------------------------")

print((x+x) == (y + y))
"""
