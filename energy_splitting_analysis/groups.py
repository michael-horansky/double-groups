
import copy

import spgrep
from improper_rotations import *


#TODO add basis functions idk


class Group():
    
    rounding_decimals = 10
    
    def __init__(self, conjugacy_classes, representations, conjugacy_class_sizes, character_table):
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
                
        self.conjugacy_class_sizes = conjugacy_class_sizes
        
        # Character table MUST be the dimension [len(representations)][len(conjugacy_classes)]
        
        self.character_table = np.array(character_table)
        
        self.conjugacy_class_elements = [] #i-th element is the list of symmetry operations in the i-th CC
        for i in range(len(self.conjugacy_classes)):
            self.conjugacy_class_elements.append([])
        
        self.group_elements = {} # a dictionary {symmetry operation : conjugacy class index}
        
    
    def add_symmetry_operation(self, group_element, conjugacy_class_index):
        self.conjugacy_class_elements[conjugacy_class_index].append(group_element)
        self.group_elements[group_element] = conjugacy_class_index
    
    def reduce_representation(self, reducible_representation):
        coefs = [0] * len(self.character_table)

        for i in range(len(self.character_table)):
            #s = 0
            for j in range(len(self.character_table[i])):
                coefs[i] += reducible_representation[j] * self.character_table[i][j] * self.conjugacy_class_sizes[j]
            coefs[i] /= sum(self.conjugacy_class_sizes)
            coefs[i] = np.round(coefs[i], decimals = Group.rounding_decimals)
            
        human_readable_output = ""
        for i in range(len(coefs)):
            if coefs[i] != 0.0:
                human_readable_output += str(int(coefs[i])) + "." + self.representations[i] + " + "

        return(coefs, human_readable_output[:-3])
            
        




# ------------------- Group generation methods ---------------------------

def name_irreps(irreps, group_operations, convention = "altmann"):
    
    # deduces the names for irreps using a specified convention
    
    final_names = []
    reorder_indices = []
    
    letter_dictionary = {
        "altmann" : [["A", "B"], "E", "T", "F", "H", "I"]
        }
    
    count_number = [[0, 0], 0, 0, 0, 0, 0]
    
    group_elements = list(group_operations.keys())
    
    # first, identify the indices of group operations that are rotations about [0, 0, 1]
    z_rot_indices = []
    for i in range(len(group_elements)):
        if np.dot(group_operations[group_elements[i]].axis, np.array([0.0, 0.0, 1.0])) == 1.0:
            z_rot_indices.append(i)
    
    
    
    def check_if_complex_conjugate(i1, i2):
        irrep1 = irreps[i1]
        irrep2 = irreps[i2]
        for i in range(len(irrep1)):
            if not np.all(irrep1[i].conjugate() == irrep2[i]):
                return(False)
        return(True)
    
    def check_if_symmetric(irrep_i):
        for index in z_rot_indices:
            if irreps[irrep_i][index][0][0] != 1.0:
                return(False)
        return(True)
    
    # first, classify irreps by their dimension
    irreps_by_dim = [] # [[dim1, [[irrep1 index], [irrep2 i, irrep2 c.c. i]...]], [dim2, [[irrep3], [irrep4]...]], ...]
    for i in range(len(irreps)):
        cur_dim = irreps[i].shape[1]
        irrep_placed = False
        biggest_smaller_dim_index = -1
        for j in range(len(irreps_by_dim)):
            if irreps_by_dim[j][0] == cur_dim:
                
                # check if there's a complex conjugate irrep present
                cc_found = False
                for k in range(len(irreps_by_dim[j][1])):
                    if len(irreps_by_dim[j][1][k]) > 1:
                        continue
                    if check_if_complex_conjugate(i, irreps_by_dim[j][1][k][0]):
                        cc_found = True
                        irrep_placed = True
                        #final_names[i] = "1" + letter_dictionary[convention][cur_dim * 2]
                        #final_names[irreps_by_dim[j][1][k][0]] = "2" + letter_dictionary[convention][cur_dim * 2]
                        #irreps_by_dim[j][1].pop(k) # removes both reps
                        irreps_by_dim[j][1][k].append(i)
                        break
                if not cc_found:
                    irreps_by_dim[j][1].append([i])
                    irrep_placed = True
                break
            elif irreps_by_dim[j][0] < cur_dim:
                biggest_smaller_dim_index = j
        if not irrep_placed:
            irreps_by_dim.insert(biggest_smaller_dim_index + 1, [cur_dim, [[i]]])
    
    for i in range(len(irreps_by_dim)):
        cur_dim = irreps_by_dim[i][0]
        for j in range(len(irreps_by_dim[i][1])):
            if len(irreps_by_dim[i][1][j]) == 1:
                cur_i = irreps_by_dim[i][1][j][0]
                reorder_indices.append(cur_i)
                if cur_dim == 1:
                    # check for antisymmetricity ALTMANN P 63
                    if check_if_symmetric(cur_i):
                        count_number[cur_dim - 1][0] += 1
                        final_names.append(letter_dictionary[convention][cur_dim-1][0] + "_" + str(count_number[cur_dim - 1][0]))
                    else:
                        count_number[cur_dim - 1][1] += 1
                        final_names.append(letter_dictionary[convention][cur_dim-1][1] + "_" + str(count_number[cur_dim - 1][1]))
                else:
                    count_number[cur_dim - 1] += 1
                    final_names.append(letter_dictionary[convention][cur_dim-1] + "_" + str(count_number[cur_dim - 1]))
            elif len(irreps_by_dim[i][1][j]) == 2:
                cur_i1 = irreps_by_dim[i][1][j][0]
                cur_i2 = irreps_by_dim[i][1][j][1]
                reorder_indices.append(cur_i1)
                reorder_indices.append(cur_i2)
                final_names.append("1" + letter_dictionary[convention][cur_dim*2-1])
                final_names.append("2" + letter_dictionary[convention][cur_dim*2-1])
    
    return(reorder_indices, final_names)
                




def group_from_multiplication_table(group_operations, multiplication_table):
    
    # group_operations is a dictionary {'label' : ImproperRotation}
    # multiplication_table is a matrix [i][j] = (element[i] . element[j]).label
    
    group_elements = list(group_operations.keys())
    
    h = len(group_elements)
    
    # First, we find the conjugacy classes
    # elements i, j are conjugate if there exists k such that multi[i][k] = multi[k][j]
    
    # For optimalization, we classify the elements by their order:
    element_order = {}
    orders_list = []
    for i in range(h):
        orders_list.append([])
    
    for i in range(h):
        j = 1
        x = group_operations[group_elements[i]]
        while(x != identity_rotation):
            j += 1
            x += group_operations[group_elements[i]]
            
            #safety break
            if j > h:
                print("ERROR: something broke in the order calculation, chief.")
                return(-1)
        element_order[group_elements[i]] = j
        orders_list[j].append(i)
    
    # Within each order list, categorize by classes
    
    def is_in_class_with(x, y):
        for k in range(h):
            if multiplication_table[x][k] == multiplication_table[k][y]:
                return(True)
        return(False)
    
    conjugacy_classes_dict = {}
    
    indices_of_representative_elements = []
    
    for order in range(h):
        headers_of_current_classes = []
        while(len(orders_list[order]) > 0):
            cur_element = orders_list[order].pop()
            element_placed = False
            for i in range(len(headers_of_current_classes)):
                if is_in_class_with(cur_element, headers_of_current_classes[i]):
                    #conjugacy_classes_dict[headers_of_current_classes[i].append(cur_element)]
                    conjugacy_classes_dict[group_elements[headers_of_current_classes[i]]].append(group_elements[cur_element])
                    element_placed = True
                    break
            if not element_placed:
                headers_of_current_classes.append(cur_element)
                conjugacy_classes_dict[group_elements[cur_element]] = [group_elements[cur_element]]
                indices_of_representative_elements.append(cur_element)
    
    #print(conjugacy_classes_dict)
    #print("------------------------------------------")
    #print(indices_of_representative_elements)
    
    
    # Second: we find the dimensions of the irreps. This is the diophantine eq sum d^2 = h, where the number of irreps = number of classes
    number_of_classes = len(conjugacy_classes_dict.keys())
    
    # the first rep has trivially d = 1
    def square_sum(l):
        res = 0
        for x in l:
            res += x * x
        return(res)
    
    target_sum = h - 1
    maximum_d = int(np.floor(np.sqrt(target_sum)))
    irrep_dimensions = [1] * (number_of_classes - 1) # monotonously increasing - check all possibilities
    while(square_sum(irrep_dimensions) != target_sum):
        index_of_increasable_irrep = len(irrep_dimensions) - 1
        while(irrep_dimensions[index_of_increasable_irrep] == maximum_d and index_of_increasable_irrep > -1):
            index_of_increasable_irrep -= 1
        if index_of_increasable_irrep == -1:
            # We have failed: no solution exists / was found
            print("ERROR: no solution to the irrep dimensionality equations was found")
            return(-1)
        
        # increase that irrep's d by 1 and set all the subsequent d's to the same value
        new_d = irrep_dimensions[index_of_increasable_irrep] + 1
        for i in range(index_of_increasable_irrep, len(irrep_dimensions)):
            irrep_dimensions[i] = new_d
    irrep_dimensions = [1] + irrep_dimensions
    print(irrep_dimensions)
    
    # at this point, we ask the user to fill out the character table themselves because difficult to automate
    
    rotation_matrices = np.zeros((len(group_elements), 3, 3))
    for i in range(len(group_elements)):
        rotation_matrices[i] = group_operations[group_elements[i]].cartesian_rep()
    
    print(rotation_matrices)
    
    irreps = spgrep.irreps.enumerate_unitary_irreps(rotation_matrices)[0]
    
    irrep_reordering, irrep_names = name_irreps(irreps, group_operations)
    
    irreps = [irreps[i] for i in irrep_reordering]
    
    for i in range(len(irreps)):
        print(irrep_names[i], irreps[i])
    
    characters = []
    for irrep in irreps:
        characters.append(spgrep.representation.get_character(irrep))
    
    conjugacy_class_names = list(conjugacy_classes_dict.keys())
    conjugacy_class_sizes = []
    for i in range(len(conjugacy_class_names)):
        conjugacy_class_sizes.append(len(conjugacy_classes_dict[conjugacy_class_names[i]]))
    
    character_table = np.zeros((len(irreps), len(irreps)))
    
    return(Group(conjugacy_class_names, irrep_names, conjugacy_class_sizes, character_table))
    
    
    #print(characters)
    
    
    # for automation, check https://github.com/gap-system/gap, https://www.gap-system.org/Overview/Capabilities/representations.html
    
    

def generate_group(generators):
    # generators is a dictionary of the form {'label' : ImproperRotation}, where the first element is 'e' : identity
    
    group_elements = list(generators.keys())
    
    group_operations = copy.deepcopy(generators)
    
    cur_h = len(group_elements)
    
    multiplication_table = []
    multiplication_table.append(group_elements.copy())
    for i in range(1, cur_h):
        multiplication_table.append([group_elements[i]] + [False]*(cur_h - 1))
    
    # We initialise the list of products which need to be checked whether they form a new group element
    products_to_check = []
    for i in range(1, cur_h):
        for j in range(1, cur_h):
            products_to_check.append([(i, j), group_operations[group_elements[i]] + group_operations[group_elements[j]]])
    
    while(len(products_to_check) > 0):
        
        cur_product = products_to_check.pop(0)
        component_i, component_j = cur_product[0]
        
        # We check if this product already exists in group elements
        exists_in_group_elements = False
        for i in range(len(group_elements)):
            if group_operations[group_elements[i]] == cur_product[1]:
                exists_in_group_elements = True
                matching_group_element_index = i
                break
        
        # If exists, we update the multiplication table and continue
        if exists_in_group_elements:
            multiplication_table[component_i][component_j] = group_elements[matching_group_element_index]
            continue
        
        #If doesn't, it forms a new group element and elongates the list
        new_label = f"{group_elements[component_i]}.{group_elements[component_j]}"
        group_operations[new_label] = cur_product[1]
        group_elements.append(new_label)
        
        # expand the multiplication table: add a column and a row
        multiplication_table[0].append(new_label)
        for i in range(1, len(multiplication_table)):
            multiplication_table[i].append(False)
        multiplication_table.append([new_label] + [False] * (len(multiplication_table[0]) - 1))
        multiplication_table[component_i][component_j] = new_label
        
        
        # Add the product "x.new" and "new.x" for all existing x in group_elements, and also "new.new"
        # Optimalization: we can ignore multiplication with e, since that is trivially already in group_elements - hence i starts at 1
        for i in range(1, len(group_elements) - 1):
            products_to_check.append([(i, len(group_elements) - 1), group_operations[group_elements[i]] + group_operations[group_elements[len(group_elements) - 1]]])
            products_to_check.append([(len(group_elements) - 1, i), group_operations[group_elements[len(group_elements) - 1]] + group_operations[group_elements[i]]])
        products_to_check.append([(len(group_elements) - 1, len(group_elements) - 1), group_operations[group_elements[len(group_elements) - 1]] + group_operations[group_elements[len(group_elements) - 1]]])
    
    # Make it into a group
    
    res_group = group_from_multiplication_table(group_operations, multiplication_table)
    
    return(group_operations, multiplication_table)


# ------------------ Example: D_4 -> D_2 -----------------

E = identity_rotation
Cz_4 = ImproperRotation([0.0, 0.0, 1.0], np.pi / 2.0, False)
Cz_2 = ImproperRotation([0.0, 0.0, 1.0], np.pi, False)
Cz_43 = ImproperRotation([0.0, 0.0, 1.0], 3.0 * np.pi / 2.0, False)

Cx_2 = ImproperRotation([1.0, 0.0, 0.0], np.pi, False)
Cy_2 = ImproperRotation([0.0, 1.0, 0.0], np.pi, False)
Cdia_1 = Cz_4 + Cy_2
Cdia_2 = Cy_2 + Cz_4

#print(Cdia_2 == Cz_4 + Cz_4 + Cz_4 + Cy_2) #LETSGOOO THIS IS TRUEEEEEEE


#print(((Cz_4 + Cy_2) + (Cz_4 + Cy_2)) == E)


operations, cayley = generate_group({"e" : E, "a" : Cz_4, "b" : Cy_2})

print(list(operations.keys()))

#print(cayley)


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
