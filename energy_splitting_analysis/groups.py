
import copy

import spgrep
from improper_rotations import *


#TODO add basis functions idk


def find_closest_matrix(matrices, m):
    i = -1
    smallest_difference = 1e9
    m_flat = np.array(m.flatten())
    for j in range(len(matrices)):
        ar_dif = (m_flat - np.array(matrices[j].flatten()))
        cur_dif = np.real(np.sum(ar_dif * np.conjugate(ar_dif)))
        if cur_dif < smallest_difference:
            smallest_difference = cur_dif
            i = j
    #print(smallest_difference)
    return(i)


def get_product_label(a, b, contraction = True):
    
    if not contraction:
        return(f"{a}.{b}")
    
    last_atom_a = a.split(".")[-1].split("^")
    first_atom_b = b.split(".")[0].split("^")
    
    
    if last_atom_a[0] == first_atom_b[0]:
        
        
        total_len_a = len(last_atom_a[0])
        total_len_b = len(first_atom_b[0])
        
        if len(last_atom_a) == 1:
            last_atom_a_power = 1
        else:
            last_atom_a_power = int(last_atom_a[1])
            total_len_a += len(last_atom_a[1]) + 1
        if len(first_atom_b) == 1:
            first_atom_b.append("")
            first_atom_b_power = 1
        else:
            first_atom_b_power = int(first_atom_b[1])
            total_len_b += len(first_atom_b[1]) + 1
        
    
        new_power = str(last_atom_a_power + first_atom_b_power)
        new_atom = last_atom_a[0] + "^" + new_power
        trim_a = a[:-total_len_a]
        trim_b = b[total_len_b:]
        
        return(trim_a + new_atom + trim_b)
    return(a + "." + b)



class Group():
    
    rounding_decimals = 10
    
    # ---------------- initializers, descriptors, destructors
    
    irrep_letter_dictionary = {
        "altmann" : [["A", "B"], "E", "T", "F", "H", "I"]
        }
    
    def __init__(self, name):
        
        self.name = name
        
        self.conjugacy_class_sizes = []
        self.irrep_dimensions = []
        
        self.generators = {}
        
        self.group_operations = {} # {"group element" : ImproperRotation / matrix} - the dictionary of the faithful representation
        self.group_elements = [] # ["group element"] - the list of names of group elements
        self.multiplication_table = [] # [i][j] = "group element" - a list of lists
        self.multiplication_dictionary = {} # {"group element" : {"group element" : "group element"}} - a dict of dicts - ONLY FOR CONVENIENCE
        
        self.z_rot_elements = [] # list of element names that are rotations around the z axis - useful for altmann naming convention
        self.indices_of_representative_elements = [] # indices in self.group_elements of the first group elements in respective conjugacy classes
        self.multiplicities = {} # {"group element" : multiplicity} - a dictionary of rotation multiplicities associated with each group element
        
        self.regular_representation = []
    
    def print_character_table(self, cc_separation = 2):
        
        def st(a, w):
            if len(a) >= w:
                return(a)
            else:
                diff = w - len(a)
                return(" " * int(np.floor(diff / 2.0)) + a + " " * int(np.ceil(diff / 2.0)))
        
        max_len_irrep = len(self.name)
        for irrep in self.irrep_names:
            if len(irrep) > max_len_irrep:
                max_len_irrep = len(irrep)
        cc_len = []
        for cc in self.conjugacy_class_names:
            cc_len.append(len(cc))
        
        printing_character_table = []
        
        for j in range(len(self.conjugacy_class_names)):
            printing_character_table.append([])
            for i in range(len(self.irrep_names)):
                if np.imag(self.character_table[j][i]) == 0.0:
                    cur_char = f"{np.real(self.character_table[j][i]):.3g}"
                    #real
                else:
                    cur_char = f"{self.character_table[j][i]:.3g}"
                printing_character_table[j].append(cur_char)
                if cc_len[i] < len(cur_char):
                    cc_len[i] = len(cur_char)
        
        header_str = st(self.name, max_len_irrep) + "|"
        for i in range(len(self.conjugacy_class_names)):
            header_str += st(self.conjugacy_class_names[i], cc_len[i] + cc_separation)
        print(header_str)
        print("-" * len(header_str))
        for i in range(len(self.irrep_names)):
            cur_str = st(self.irrep_names[i], max_len_irrep) + "|"
            for j in range(len(self.conjugacy_class_names)):
                cur_str += st(printing_character_table[i][j], cc_len[j] + cc_separation)
            print(cur_str)
            
    
    # --------------- property management methods
    
    
    def set_conjugacy_classes(self, arg1, arg2 = -1):
        
        # Full specification:    dict {"name" : [elements]}
        # Minimal specification: (sizes, names). If no names, names become labelled by sizes and index
        
        if type(arg1) == dict:
            self.conjugacy_classes = arg1.copy()
            self.conjugacy_class_names = list(self.conjugacy_classes.keys())
            self.conjugacy_class_sizes = []
            for i in range(len(self.conjugacy_class_names)):
                self.conjugacy_class_sizes.append(len(self.conjugacy_classes[self.conjugacy_class_names[i]]))
        
        else:
            self.conjugacy_classes = {}
            self.conjugacy_class_sizes = arg1.copy()
            if arg2 != -1:
                self.conjugacy_class_names = arg2.copy()
            else:
                self.conjugacy_class_names = []
                for i in range(len(self.conjugacy_class_sizes)):
                    self.conjugacy_class_names.append(str(self.conjugacy_class_sizes[i]) + "A_" + str(i + 1))
        
    
    def set_irreducible_representations(self, arg1, arg2 = -1):
        
        # Full specification: dict {"name" : ndarray[matrix 1, matrix 2, ...]}
        # Minimal specification: ([dimensions], [names]). if no names, names becomes labelled minimally according to altmann
        
        if type(arg1) == dict:
            self.irreps = arg1.copy()
            self.irrep_names = list(self.irreps.keys())
            self.irrep_dimensions = []
            for i in range(len(self.irrep_names)):
                self.irrep_dimensions.append(self.irreps[self.irrep_names[i]].shape[1])
        
        elif type(arg1[0]) == np.ndarray or type(arg1[0]) == list:
            # unnamed but full irrep knowledge
            self.name_and_set_irreps(arg1)
        
        else:
            self.irreps = {}
            self.irrep_dimensions = arg1.copy()
            if arg2 != -1:
                self.irrep_names = arg2.copy()
            else:
                self.irrep_names = []
                count_number = [0, 0, 0, 0, 0, 0]
                
                for i in range(len(self.irrep_dimensions)):
                    if self.irrep_dimensions[i] == 1:
                        letter = f"{Group.irrep_letter_dictionary['altmann'][0][0]}/{Group.irrep_letter_dictionary['altmann'][0][1]}"
                    else:
                        letter = Group.irrep_letter_dictionary['altmann'][irrep_dimensions[i] - 1]
                    self.irrep_names.append(letter + "_" + str(i + 1))
        
    def set_character_table(self, character_table):
        
        # Character table MUST be the dimension [len(representations)][len(conjugacy_classes)]
        # either a dict of dicts {"irrep name" : {"conjugacy class name" : character}} -> we transform to list of lists that is well ordered
        # or a list of lists [irrep index][cc index] = character (THIS ASSUMES THE IRREPS ARE ORDERED IN THE SAME WAY AS THE ASSIGNED IRREPS; SAME FOR CCs)
        
        # If irrep dimensions are an empty list - no irreps specified, their order is determined by the character table
        # We always assume that conjugacy classes are initialized, with the first one being the identity class
        if self.irrep_dimensions == []:
            if type(character_table) == dict:
                self.irrep_names = list(character_table.keys())
                for i in range(len(self.irrep_names)):
                    self.irrep_dimensions.append(int(character_table[self.irrep_names[i]][self.conjugacy_class_names[0]]))
            else:
                for i in range(len(character_table)):
                    self.irrep_dimensions.append(int(character_table[i][0]))
                self.irrep_names = []
                count_number = [0, 0, 0, 0, 0, 0]
                for i in range(len(self.irrep_dimensions)):
                    if self.irrep_dimensions[i] == 1:
                        letter = f"{Group.irrep_letter_dictionary['altmann'][0][0]}/{Group.irrep_letter_dictionary['altmann'][0][1]}"
                    else:
                        letter = Group.irrep_letter_dictionary['altmann'][irrep_dimensions[i] - 1]
                    self.irrep_names.append(letter + "_" + str(i + 1))
        
        
        self.character_table = np.zeros((len(self.irrep_names), len(self.conjugacy_class_names)))
        
        if type(character_table) == dict:
            for i in range(len(self.irrep_names)):
                for j in range(len(self.conjugacy_class_names)):
                    self.character_table[i][j] = character_table[self.irrep_names[i]][self.conjugacy_class_names[j]]
        
        else:
            for i in range(len(self.irrep_names)):
                for j in range(len(self.conjugacy_class_names)):
                    self.character_table[i][j] = character_table[i][j]
    
    def name_and_set_irreps(self, irreps, convention = "altmann"):
    
        # deduces the names for irreps using a specified convention and saves them as a dictionary
        # this function can be called INSTEAD of self.set_irreducible_representations()
        
        final_names = []
        reorder_indices = []
        
        count_number = [[0, 0], 0, 0, 0, 0, 0]
        
        characters = []
        for irrep in irreps:
            characters.append(spgrep.representation.get_character(irrep))
        
        #group_elements = list(group_operations.keys())
        """
        # first, identify the indices of group operations that are rotations about [0, 0, 1]
        z_rot_indices = []
        
        # This is either SO'3 or SU 2
        if faithful_rep.shape[1] == 3:
            # SO'3
            for i in range(len(group_elements)):
                if np.all(np.matmul(faithful_rep[i], np.array([0.0, 0.0, 1.0])) == np.array([0.0, 0.0, 1.0])):
                    z_rot_indices.append(i)
        elif faithful_rep.shape[1] == 2:
            # SU 2
            for i in range(len(group_elements)):
                if faithful_rep[i][0][0] == faithful_rep[i][1][1] and faithful_rep[i][0][1] == faithful_rep[i][1][0]:
                    z_rot_indices.append(i)"""
        """for i in range(len(group_elements)):
            if np.dot(group_operations[group_elements[i]].axis, np.array([0.0, 0.0, 1.0])) == 1.0:
                z_rot_indices.append(i)"""
        
        if self.z_rot_elements == []:
            # we don't know which elements are rotations about the z axis - for 1D reps we use "A/B"
            z_rot_indices = []
        else:
            z_rot_indices = []
            for z_rot_element in self.z_rot_elements:
                z_rot_indices.append(self.group_elements.index(z_rot_element))
        
        def check_if_complex_conjugate(i1, i2):
            # deducing from the T group character table in Altmann, this refers to a relation between the _characters_ of the irreps
            
            return(np.all(np.round(characters[i1].conjugate(), 5) == np.round(characters[i2], 5)))
            """
            irrep1 = irreps[i1]
            irrep2 = irreps[i2]
            for i in range(len(irrep1)):
                if not np.all(irrep1[i].conjugate() == irrep2[i]):
                    return(False)
            return(True)"""
        
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
                        if self.z_rot_elements == []:
                            count_number[cur_dim - 1][0] += 1
                            count_number[cur_dim - 1][1] += 1
                            cur_letter = f"{Group.irrep_letter_dictionary[convention][cur_dim-1][0]}/{Group.irrep_letter_dictionary[convention][cur_dim-1][1]}"
                            final_names.append(cur_letter + "_" + str(count_number[cur_dim - 1][0]))
                        else:
                            if check_if_symmetric(cur_i):
                                count_number[cur_dim - 1][0] += 1
                                final_names.append(Group.irrep_letter_dictionary[convention][cur_dim-1][0] + "_" + str(count_number[cur_dim - 1][0]))
                            else:
                                count_number[cur_dim - 1][1] += 1
                                final_names.append(Group.irrep_letter_dictionary[convention][cur_dim-1][1] + "_" + str(count_number[cur_dim - 1][1]))
                    else:
                        count_number[cur_dim - 1] += 1
                        final_names.append(Group.irrep_letter_dictionary[convention][cur_dim-1] + "_" + str(count_number[cur_dim - 1]))
                elif len(irreps_by_dim[i][1][j]) == 2:
                    cur_i1 = irreps_by_dim[i][1][j][0]
                    cur_i2 = irreps_by_dim[i][1][j][1]
                    reorder_indices.append(cur_i1)
                    reorder_indices.append(cur_i2)
                    final_names.append("1" + Group.irrep_letter_dictionary[convention][cur_dim*2-1])
                    final_names.append("2" + Group.irrep_letter_dictionary[convention][cur_dim*2-1])
        
        # reorder irreps, since dict remembers insertion order
        
        irrep_dict = {}
        for i in range(len(reorder_indices)):
            j = reorder_indices[i]
            irrep_dict[final_names[j]] = irreps[j]
        self.set_irreducible_representations(irrep_dict)
        
        
        
    # ----------------- group generation methods
    
    def generate_multiplication_table(self, generators):
        
        # REQUIREMENTS: none
        # SPECIFIED PARAMETERS: self.generators, self.group_elements, self.group_operations, self.multiplication_table, self.z_rot_elements
        
        # generators is a dictionary of the form {'label' : ImproperRotation}, where the first element is 'e' : identity
        self.generators = generators
        
        self.group_elements = list(generators.keys())
        
        self.group_operations = copy.deepcopy(generators)
        
        cur_h = len(self.group_elements)
        
        self.multiplication_table = []
        self.multiplication_table.append(self.group_elements.copy())
        for i in range(1, cur_h):
            self.multiplication_table.append([self.group_elements[i]] + [False]*(cur_h - 1))
        
        # We initialise the list of products which need to be checked whether they form a new group element
        products_to_check = []
        for i in range(1, cur_h):
            for j in range(1, cur_h):
                products_to_check.append([(i, j), self.group_operations[self.group_elements[i]] + self.group_operations[self.group_elements[j]]])
        
        while(len(products_to_check) > 0):
            
            #print(self.group_elements)
            
            cur_product = products_to_check.pop(0)
            component_i, component_j = cur_product[0]
            
            # We check if this product already exists in group elements
            exists_in_group_elements = False
            for i in range(len(self.group_elements)):
                if self.group_operations[self.group_elements[i]] == cur_product[1]:
                    exists_in_group_elements = True
                    matching_group_element_index = i
                    break
            
            # If exists, we update the multiplication table and continue
            if exists_in_group_elements:
                self.multiplication_table[component_i][component_j] = self.group_elements[matching_group_element_index]
                continue
            
            #If doesn't, it forms a new group element and elongates the list
            new_label = get_product_label(self.group_elements[component_i], self.group_elements[component_j])
            self.group_operations[new_label] = cur_product[1]
            self.group_elements.append(new_label)
            
            # expand the multiplication table: add a column and a row
            self.multiplication_table[0].append(new_label)
            for i in range(1, len(self.multiplication_table)):
                self.multiplication_table[i].append(False)
            self.multiplication_table.append([new_label] + [False] * (len(self.multiplication_table[0]) - 1))
            self.multiplication_table[component_i][component_j] = new_label
            
            
            # Add the product "x.new" and "new.x" for all existing x in group_elements, and also "new.new"
            # Optimalization: we can ignore multiplication with e, since that is trivially already in group_elements - hence i starts at 1
            for i in range(1, len(self.group_elements) - 1):
                products_to_check.append([(i, len(self.group_elements) - 1), self.group_operations[self.group_elements[i]] + self.group_operations[self.group_elements[len(self.group_elements) - 1]]])
                products_to_check.append([(len(self.group_elements) - 1, i), self.group_operations[self.group_elements[len(self.group_elements) - 1]] + self.group_operations[self.group_elements[i]]])
            products_to_check.append([(len(self.group_elements) - 1, len(self.group_elements) - 1), self.group_operations[self.group_elements[len(self.group_elements) - 1]] + self.group_operations[self.group_elements[len(self.group_elements) - 1]]])
        
        # now, we check which elements are rotations around the z axis, and commit multiplicities to memory:
        self.z_rot_elements = []
        self.multiplicities = {}
        for group_element in self.group_elements:
            if np.all(np.matmul(self.group_operations[group_element].cartesian_rep(), np.array([0.0, 0.0, 1.0])) == np.array([0.0, 0.0, 1.0])):
                self.z_rot_elements.append(group_element)
            self.multiplicities[group_element] = self.group_operations[group_element].multiplicity
        
        self.order = len(self.group_elements)
        

    def get_double_group(self):
        
        # REQUIREMENTS: self.group_elements, self.group_operations, self.z_rot_elements
        # SPECIFIED PARAMETERS: self.group_elements, self.group_operations, self.z_rot_elements, self.multiplication_table
        
        # This is the group doubling method. A few notes:
        #   -This method overwrites self.group_elements, self.group_operations, and self.multiplication_table, so dependent methods have to be re-called
        #   -The resulting double group's faithful rep is no longer expressable as improper rotations, but rather as SU(2) matrices + inversions.
        
        element_matrices = [0] * self.order * 2
        
        for i in range(self.order):
            element_matrices[i]     = self.group_operations[self.group_elements[i]].SU2_rep()
            element_matrices[i + self.order] = - self.group_operations[self.group_elements[i]].SU2_rep()
            if i == 0:
                new_label = "R"
            else:
                new_label = "R" + self.group_elements[i]
            self.group_elements.append(new_label)
            if self.group_elements[i] in self.z_rot_elements:
                self.z_rot_elements.append(new_label)
            
        
        # now we generate and populate the new multiplication table
        
        new_mt = []
        for i in range(2 * self.order):
            new_mt.append([])
            for j in range(2 * self.order):
                product_m = np.matmul(element_matrices[i], element_matrices[j])
                product_index = find_closest_matrix(element_matrices, product_m)
                new_mt[i].append(self.group_elements[product_index])
        
        self.multiplication_table = new_mt
        # update the group operations and multiplicities
        # the multiplicity of a time reversal rotation is unconstrained - angle is bigger than 2pi
        for i in range(len(self.group_elements)):
            self.group_operations[self.group_elements[i]] = element_matrices[i]
        for i in range(self.order):
            p = self.multiplicities[self.group_elements[i]][0]
            q = self.multiplicities[self.group_elements[i]][1]
            self.multiplicities[self.group_elements[i + self.order]] = [p + q, q]
        
        # Now we update the group order
        self.order = len(self.group_elements)
            

    
    
    def conjugacy_classes_from_multiplication_table(self):
        # REQUIREMENTS: self.group_elements, self.multiplication_table
        # SPECIFIED PARAMETERS: self.conjugacy_classes, self.conjugacy_class_names, self.conjugacy_class_sizes
        # this function can be called INSTEAD of self.set_conjugacy_classes()
        
        self.order = len(self.group_elements)
        # elements i, j are conjugate if there exists k such that multi[i][k] = multi[k][j]
        
        # For optimalization, we classify the elements by their order:
        element_order = {}
        orders_list = []
        for i in range(self.order):
            orders_list.append([])
        
        for i in range(self.order):
            j = 1
            #x = group_operations[self.group_elements[i]]
            x = i
            while(x != 0):
                j += 1
                #x += group_operations[self.group_elements[i]]
                x = self.group_elements.index(self.multiplication_table[x][i])
                
                #safety break
                if j > self.order:
                    print("ERROR: something broke in the order calculation, chief.")
                    return(-1)
            element_order[self.group_elements[i]] = j
            orders_list[j].append(i)
        
        # Within each order list, categorize by classes
        
        def is_in_class_with(x, y):
            for k in range(self.order):
                if self.multiplication_table[x][k] == self.multiplication_table[k][y]:
                    return(True)
            return(False)
        
        conjugacy_classes_dict = {}
        
        # TODO: conjugacy_class_names, add x+Rx naming convention
        
        self.indices_of_representative_elements = []
        
        for order in range(self.order):
            headers_of_current_classes = []
            while(len(orders_list[order]) > 0):
                cur_element = orders_list[order].pop()
                element_placed = False
                for i in range(len(headers_of_current_classes)):
                    if is_in_class_with(cur_element, headers_of_current_classes[i]):
                        #conjugacy_classes_dict[headers_of_current_classes[i].append(cur_element)]
                        conjugacy_classes_dict[self.group_elements[headers_of_current_classes[i]]].append(self.group_elements[cur_element])
                        element_placed = True
                        break
                if not element_placed:
                    headers_of_current_classes.append(cur_element)
                    conjugacy_classes_dict[self.group_elements[cur_element]] = [self.group_elements[cur_element]]
                    self.indices_of_representative_elements.append(cur_element)
        
        # TODO rename dictionary keys here (use d['new'] = d.pop('old'))
        # Add conjugacy class sizes and x + Rx to relevant ccs
        
        self.set_conjugacy_classes(conjugacy_classes_dict)
    
    def get_character_table(self, tolerance_decimals = 5):
        
        # REQUIREMENTS: self.irreps, self.irrep_names, self.group_elements, self.conjugacy_class_names, self.indices_of_representative_elements
        # SPECIFIED PARAMETERS: self.character_table
        # this function can be called INSTEAD of self.set_character_table()
        # to have the indices of representative elements, this function should only be called after self.conjugacy_classes_from_multiplication_table
        
        characters = []
        for irrep in self.irrep_names:
            characters.append(spgrep.representation.get_character(self.irreps[irrep]))
        
        #conjugacy_class_sizes = []
        #for i in range(len(conjugacy_class_names)):
        #    conjugacy_class_sizes.append(len(conjugacy_classes_dict[conjugacy_class_names[i]]))
        
        self.character_table = np.zeros((len(self.irreps), len(self.irreps)), dtype=np.complex_)
        
        for i_irrep in range(len(self.irreps)):
            for i_cc in range(len(self.conjugacy_class_names)):
                self.character_table[i_irrep][i_cc] = np.round(characters[i_irrep][self.indices_of_representative_elements[i_cc]], decimals = tolerance_decimals)
    
    def get_regular_representation(self):
        
        # REQUIREMENTS: self.group_elements, self.multiplication_table
        # SPECIFIED PARAMETERS: self.regular_representation
        
        # assume the neutral element is first in group_elements
        self.order = len(self.group_elements)
        group_element_indices = {}
        for i in range(self.order):
            group_element_indices[self.group_elements[i]] = i
        
        # First, we rearrange the mult_table so that E is along the main diagonal
        mt = np.array(self.multiplication_table)
        # we do this by swapping rows
        for i in range(1, self.order):
            # we want to find the row which has E as its i-th element
            if mt[i][i] == self.group_elements[0]:
                continue
            for j in range(i + 1, self.order + 1):
                if j == self.order:
                    print("ERROR: multiplication table isn't proper")
                    #continue
                    return(-1)
                if mt[j][i] == self.group_elements[0]:
                    # swap i and j rows in mt
                    mt[[i, j]] = mt[[j, i]]
                    break
        
        # By swapping rows, the first index of mt no longer corresponds to its group_element item, but rather its inverse
        
        self.regular_representation = np.zeros((self.order, self.order, self.order))
        for i in range(self.order):
            for j in range(self.order):
                self.regular_representation[group_element_indices[mt[i][j]]][i][j] = 1.0
    
    
    def initialize_from_multiplication_table(self):
        
        # REQUIREMENTS: self.group_elements, self.multiplication_table
        # SPECIFIED PARAMETERS: everything
        
        self.order = len(self.group_elements)
        
        # First, we find the conjugacy classes
        self.conjugacy_classes_from_multiplication_table()
        
        # Here we find the regular representation
        self.get_regular_representation()
        #for i in range(len(group_elements)):
        #    print(group_elements[i], regular_rep[i])
        
        # We obtain the irreps
        irreps = spgrep.irreps.enumerate_unitary_irreps_from_regular_representation(self.regular_representation)
        for i in range(len(irreps)):
            irreps[i] = np.round(irreps[i], decimals = 12)
        
        # We name and classify irreps
        self.name_and_set_irreps(irreps)
        
        
        #irrep_reordering, irrep_names = name_irreps(irreps, group_elements, rotation_matrices)
        #irreps = [irreps[i] for i in irrep_reordering]
        
        # We set the character table
        self.get_character_table()
        
        # for automation, check https://github.com/gap-system/gap, https://www.gap-system.org/Overview/Capabilities/representations.html
    
    def generate_group(self, generators):
        
        # A wrapper for self.generate_multiplication_table() and self.initialize_from_multiplication_table()
        self.generate_multiplication_table(generators)
        self.initialize_from_multiplication_table()
    
    
    def generate_double_group(self, generators):
        
        # A wrapper for self.generate_multiplication_table(), self.get_double_group(), and self.initialize_from_multiplication_table()
        self.generate_multiplication_table(generators)
        self.get_double_group()
        self.initialize_from_multiplication_table()
        
        
    
    
    
    def add_symmetry_operation(self, group_element, conjugacy_class_index):
        self.conjugacy_class_elements[conjugacy_class_index].append(group_element)
        self.group_elements[group_element] = conjugacy_class_index
    
    def reduce_representation(self, reducible_representation):
        
        # This is either a dict {"conjugacy class name" : character} or a list [cc_i] = character
        
        coefs = [0] * len(self.character_table)
        
        if type(reducible_representation) == dict:
            for i in range(len(self.character_table)):
                for j in range(len(self.character_table[i])):
                    coefs[i] += reducible_representation[self.conjugacy_class_names[j]] * self.character_table[i][j] * self.conjugacy_class_sizes[j]
                coefs[i] /= self.order
                coefs[i] = np.round(coefs[i], decimals = Group.rounding_decimals)
                
                if np.imag(coefs[i]) != 0.0:
                    print("CAREFUL! The input rep has imaginary coefficients in its reduction; this has been omitted, but requires manual checking!!!")
                coefs[i] = np.real(coefs[i])
                
        else:
            for i in range(len(self.character_table)):
                #s = 0
                for j in range(len(self.character_table[i])):
                    coefs[i] += reducible_representation[j] * self.character_table[i][j] * self.conjugacy_class_sizes[j]
                coefs[i] /= sum(self.conjugacy_class_sizes)
                coefs[i] = np.round(coefs[i], decimals = Group.rounding_decimals)
                
                if np.imag(coefs[i]) != 0.0:
                    print("CAREFUL! The input rep has imaginary coefficients in its reduction; this has been omitted, but requires manual checking!!!")
                coefs[i] = np.real(coefs[i])
                
        human_readable_output = ""
        for i in range(len(coefs)):
            if coefs[i] != 0.0:
                if coefs[i] == 1:
                    human_readable_output += self.irrep_names[i] + " + "
                else:
                    human_readable_output += str(int(coefs[i])) + "." + self.irrep_names[i] + " + "

        return(coefs, human_readable_output[:-3])
    
    def angular_representation(self, j, symmetry = "g"):
        # for a given j value, this creates the reducible angular representation (as a subgroup of the full rotation group)
        # if symmetry = "g" (gerade), then characters dont flip sign for inversions; for symmetry = "u" (ungerade), they do
        
        rep = {}
        for i in range(len(self.conjugacy_class_names)):
            cur_mult = self.multiplicities[self.group_elements[self.indices_of_representative_elements[i]]]
            cur_angle = 2.0 * np.pi * cur_mult[0] / cur_mult[1]
            if cur_angle == 0:
                # identity class - this gives (j + 1/2) / (1/2) = 2j + 1
                rep[self.conjugacy_class_names[i]] = 2 * j + 1
            elif cur_angle == 2.0 * np.pi:
                # time reversal class - for integer j, this is 2j+1; for half-integer j, this is -(2j+1)
                if (2 * j) % 2 == 1:
                    rep[self.conjugacy_class_names[i]] = - (2 * j + 1)
                else:
                    rep[self.conjugacy_class_names[i]] = 2 * j + 1
            else:
                rep[self.conjugacy_class_names[i]] = np.sin((j + 0.5) * cur_angle) / np.sin(0.5 * cur_angle)
        return(rep)
        


# ------------------- Group generation methods ---------------------------
                

def find_irrep_dimensions(conjugacy_classes_dict):
    # We find the dimensions of the irreps. This is the diophantine eq sum d^2 = h, where the number of irreps = number of classes
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
    return(irrep_dimensions)




def group_from_group_operations(group_operations, multiplication_table):
    
    # group_operations is a dictionary {'label' : ImproperRotation}
    # multiplication_table is a matrix [i][j] = (element[i] . element[j]).label
    
    group_elements = list(group_operations.keys())
    
    h = len(group_elements)
    
    # First, we find the conjugacy classes
    conjugacy_classes_dict, indices_of_representative_elements = conjugacy_classes_from_multiplication_table(group_elements, multiplication_table)
    
    rotation_matrices = np.zeros((len(group_elements), 3, 3))
    for i in range(len(group_elements)):
        rotation_matrices[i] = group_operations[group_elements[i]].cartesian_rep()
    
    irreps = spgrep.irreps.enumerate_unitary_irreps(rotation_matrices)[0]
    
    irrep_reordering, irrep_names = name_irreps(irreps, group_elements, rotation_matrices)
    irreps = [irreps[i] for i in irrep_reordering]
    
    character_table = get_character_table(irreps, conjugacy_classes_dict, indices_of_representative_elements)
    
    conjugacy_class_names = list(conjugacy_classes_dict.keys())
    conjugacy_class_sizes = []
    for i in range(len(conjugacy_class_names)):
        conjugacy_class_sizes.append(len(conjugacy_classes_dict[conjugacy_class_names[i]]))
        conjugacy_class_names[i] = str(conjugacy_class_sizes[i]) + conjugacy_class_names[i]
    
    return(Group(conjugacy_class_names, irrep_names, conjugacy_class_sizes, character_table))
    # for automation, check https://github.com/gap-system/gap, https://www.gap-system.org/Overview/Capabilities/representations.html
    







# ------------------ Example: D_4 -> D_2 -----------------

E = identity_rotation
"""
Cz_4 = ImproperRotation([0.0, 0.0, 1.0], np.pi / 2.0, False)
Cz_2 = ImproperRotation([0.0, 0.0, 1.0], np.pi, False)
Cz_3 = ImproperRotation([0.0, 0.0, 1.0], 2.0 * np.pi / 3.0, False)
Cz_43 = ImproperRotation([0.0, 0.0, 1.0], 3.0 * np.pi / 2.0, False)
Cz_5 = ImproperRotation([0.0, 0.0, 1.0], 2.0 * np.pi / 5.0, False)

Cx_2 = ImproperRotation([1.0, 0.0, 0.0], np.pi, False)
Cy_2 = ImproperRotation([0.0, 1.0, 0.0], np.pi, False)
"""
Cz_4 = ImproperRotation([0.0, 0.0, 1.0], [1, 4], False)
Cz_2 = ImproperRotation([0.0, 0.0, 1.0], [1, 2], False)
Cz_3 = ImproperRotation([0.0, 0.0, 1.0], [1, 3], False)
Cz_43 = ImproperRotation([0.0, 0.0, 1.0], [3, 4], False)
Cz_5 = ImproperRotation([0.0, 0.0, 1.0], [1, 5], False)

Cx_2 = ImproperRotation([1.0, 0.0, 0.0], [1, 2], False)
Cy_2 = ImproperRotation([0.0, 1.0, 0.0], [1, 2], False)
Cdia_1 = Cz_4 + Cy_2
Cdia_2 = Cy_2 + Cz_4


a = np.matmul(Cz_4.SU2_rep(), Cz_4.SU2_rep())
b = -Cz_2.SU2_rep()

#print(find_closest_matrix([Cz_4.SU2_rep(), a], b))


"""print("Cz_2 =", Cz_2.axis, Cz_2.multiplicity)
print(u.axis, u.multiplicity)
print(v.axis, v.multiplicity)

print(inverse_angle([1, 4]))

print(u == v)"""

"""
print(Cx_2.SU2_rep())
print(Cz_5.SU2_rep())
print((Cx_2 + Cz_5).SU2_rep())
print(np.matmul(Cx_2.SU2_rep(), Cz_5.SU2_rep()))
print(np.all(np.matmul(Cx_2.SU2_rep(), Cz_5.SU2_rep()) == (Cx_2 + Cz_5).SU2_rep())) #WAHOOO
print("-------------------------------")
print(np.matmul(Cx_2.SU2_rep(), Cx_2.SU2_rep()))
print(np.all(np.matmul(Cx_2.SU2_rep(), Cx_2.SU2_rep()) == (Cx_2 + Cx_2).SU2_rep())) #WAHOOO
"""

#print(Cdia_2 == Cz_4 + Cz_4 + Cz_4 + Cy_2) #LETSGOOO THIS IS TRUEEEEEEE


#print(((Cz_4 + Cy_2) + (Cz_4 + Cy_2)) == E)


#operations, cayley = generate_group({"e" : E, "a" : Cz_4, "b" : Cy_2})

"""
D4_group = Group("D4")
D4_group.generate_double_group({"e" : E, "a" : Cz_4, "b" : Cy_2})
#print(D4_group.conjugacy_classes)
print(D4_group.irrep_names)
print(D4_group.character_table)


D6_group = generate_double_group({"E" : E, "Cz_3" : Cz_3, "Cz_2" : Cz_2, "C'_2" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], False)})
print(D6_group.conjugacy_classes)
print(D6_group.representations)
print(D6_group.character_table)"""

T_group = Group("T")
T_group.generate_double_group({"E" : E, "Cz_2" : Cz_2, "Cy_2" : Cy_2, "C'_3" : ImproperRotation([1.0, 1.0, 1.0], [1, 3], False)})
#print(T_group.conjugacy_class_names)
#print(T_group.irrep_names)
#print(T_group.character_table)
T_group.print_character_table()

a = T_group.angular_representation(5/2)
print(T_group.reduce_representation(a))


#TODO - direct product of groups
#TODO - symmetry breakage

"""
C3v_group = generate_group({"E" : E, "Cz_3" : Cz_3, "s" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
print(C3v_group.conjugacy_classes)
print(C3v_group.representations)
print(C3v_group.character_table)"""


#print(list(operations.keys()))

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
