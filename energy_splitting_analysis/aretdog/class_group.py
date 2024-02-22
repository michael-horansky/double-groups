# ---------------------------------------------------------
#
#   The leading src file of the AReTDoG package
#   AReTDoG: Algorithmic Representation Theory of Double Groups
#   Import this file to access the full scope of AReTDoG
#
#   This file provides:
#       class Group()
# ---------------------------------------------------------
#   Created: 3 November 2023
#   Author: Michal Horansky (michal.horansky20@imperial.ac.uk)
#
# ---------------------------------------------------------


import copy

import spgrep
from .class_improper_rotation import *
from .class_representation import *


#TODO add basis functions idk





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
        self.irrep_dimensions = {}
        
        self.generators = {}
        
        self.group_operations = {} # {"group element" : ImproperRotation / matrix} - the dictionary of the faithful representation
        self.group_elements = [] # ["group element"] - the list of names of group elements
        self.multiplication_table = [] # [i][j] = "group element" - a list of lists
        self.multiplication_dictionary = {} # {"group element" : {"group element" : "group element"}} - a dict of dicts - ONLY FOR CONVENIENCE
        self.character_table = [] # [irrep index][cc index] = character
        self.irrep_characters = {} # {"irrep name" : Representation}
        
        self.z_rot_elements = [] # list of element names that are rotations around the z axis - useful for altmann naming convention
        self.indices_of_representative_elements = {} # {"cc name" : index in self.group_elements of the first group elements in respective conjugacy class}
        self.conjugacy_class_time_reversal = {} # {"conjugacy class name" = is this class associated with time reversal?}
        self.example_properness_conjugacy_class_pair = ["E"] # pair of two conjugacy classes that differ only in properness
        self.element_spatial_properties = {} # {"group element" : [axis, multiplicity, time reversal, inversion]} - a dictionary of properties associated with each group element
        self.element_conjugacy_classes = {} # dict {"group element" : conjugacy class name}
        
        self.conjugacy_class_associated_angle = {} # {"conjugacy class name" : angle such that char(E_1/2) = 2 cos (angle/2) times time reversal}
        
        self.cartesian_basis = {}
        
        self.regular_representation = []
        
        # group operation properties
        self.subgroups = [] # a list of subgroup names
        self.subgroup_element_relations = {} # dict {"subgroup name" : [instance of group, {"subgroup el." : "group el."}]}
        self.subgroup_conjugacy_relations = {} # dict {"subgroup name" : [instance of group, {"subgroup cc." : "group cc."}]}
        
        # self descriptors
        self.is_proper = True
        self.is_double = False
        self.is_inversion_symmetry_determinable = False
        self.identity_irrep = "" # name
        self.inversion_irrep = "" # this is the label of the inversion irrep (char(proper) = 1, char(rotoinversion) = -1). For proper groups this is the identity rep.
        self.complex_conjugate_irreps = {} # {"1gamma" : "2gamma"}
    
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
        
        basis_column_name = "basis"
        basis_function_labels = []
        for irrep in self.irrep_names:
            if irrep in self.cartesian_basis.keys():
                basis_function_labels.append("(" + "); (".join(self.cartesian_basis[irrep]) + ")")
            else:
                basis_function_labels.append("")
        cc_len.append(len(basis_column_name))
        
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
            printing_character_table[j].append(basis_function_labels[j])
            if cc_len[-1] < len(basis_function_labels[j]):
                cc_len[-1] = len(basis_function_labels[j])
        
        header_str = st(self.name, max_len_irrep + cc_separation) + "|"
        for i in range(len(self.conjugacy_class_names)):
            header_str += st(self.conjugacy_class_names[i], cc_len[i] + cc_separation)
        header_str += "|" + st(basis_column_name, cc_len[-1] + cc_separation)
        print(header_str)
        print("-" * len(header_str))
        for i in range(len(self.irrep_names)):
            cur_str = st(self.irrep_names[i], max_len_irrep + cc_separation) + "|"
            for j in range(len(self.conjugacy_class_names)):
                cur_str += st(printing_character_table[i][j], cc_len[j] + cc_separation)
            cur_str += "|" + st(printing_character_table[i][-1], cc_len[-1] + cc_separation)
            print(cur_str)
        print("-" * len(header_str))
            
    
    # --------------- property management methods
    
    def set_element_conjugacy_classes(self):
        self.element_conjugacy_classes = {}
        for cc in self.conjugacy_class_names:
            for element in self.conjugacy_classes[cc]:
                self.element_conjugacy_classes[element] = cc
    
    def set_conjugacy_classes(self, arg1, arg2 = -1):
        
        # Full specification:    dict {"name" : [elements]}
        # Minimal specification: (sizes, names). If no names, names become labelled by sizes and index
        
        if type(arg1) == dict:
            self.conjugacy_classes = arg1.copy()
            self.conjugacy_class_names = list(self.conjugacy_classes.keys())
            self.conjugacy_class_sizes = []
            for i in range(len(self.conjugacy_class_names)):
                self.conjugacy_class_sizes.append(len(self.conjugacy_classes[self.conjugacy_class_names[i]]))
            self.set_element_conjugacy_classes()
        
        else:
            self.conjugacy_classes = {}
            self.conjugacy_class_sizes = arg1.copy()
            if arg2 != -1:
                self.conjugacy_class_names = arg2.copy()
            else:
                self.conjugacy_class_names = []
                for i in range(len(self.conjugacy_class_sizes)):
                    self.conjugacy_class_names.append(str(self.conjugacy_class_sizes[i]) + "A_" + str(i + 1))
        
        self.indices_of_representative_elements = {}
        for cc in self.conjugacy_class_names:
            self.indices_of_representative_elements[cc] = self.group_elements.index(self.conjugacy_classes[cc][0])
        self.set_element_conjugacy_classes()
    
    
    def classify_conjugacy_classes(self, conjugacy_class_proposal):
        # either a dict (with provisional names) or a list
        # this method a) names all conjugacy classes in a sane way, b) determines their time-reversality, and c) changes self.indices_of_representative_elements
        
        def is_contracted(cc):
            # a 'contracted' cc contains Rg for each g contained therein
            # this is uniquely determined by whether this is true for any single element - we check for the first one
            return("R" + cc[0] in cc or cc[0][1:] in cc) # a dirty check! depends on our naming convention
        
        def find_reversed_class(cc, other_classes):
            # for a non-contracted cc, finds a cc which contains Rg for all g in the original cc
            # this is uniquely determined by any single element
            for j in range(len(other_classes)):
                if "R" + cc[0] in other_classes[j] or cc[0][1:] in other_classes[j]:
                    return(j)
        
        def get_nondescript_desc(notation, marker = "'"):
            if notation in number_of_unclassifiable_multiplicities:
                number_of_unclassifiable_multiplicities[notation] += 1
                return(marker * number_of_unclassifiable_multiplicities[notation])
            else:
                number_of_unclassifiable_multiplicities[notation] = 1
                return(marker)
                
        
        sanitized_conjugacy_classes = {}
        
        number_of_unclassifiable_multiplicities = {}
        if type(conjugacy_class_proposal) == dict:
            conjugacy_class_proposal = list(conjugacy_class_proposal.values())
            
        self.is_inversion_symmetry_determinable = False # only true if there exists an inversion in a class by itself
        
        skip_classes = [] # if we take care of a reversed class, we can skip it when encountering its index
        time_reversal_pairs = {}
        for cc_i in range(len(conjugacy_class_proposal)):
            if cc_i in skip_classes:
                continue
            cc = conjugacy_class_proposal[cc_i]
            # First, the size
            size = len(cc)
            rep_spatial = self.element_spatial_properties[cc[0]]
            
            # If cc of size 1, we can label it by its only element (takes care of identities, time reversals, inversions, and mirrors)
            if size == 1:
                if not rep_spatial[2] and rep_spatial[3]:
                    # elementary improper operation, which always exists in a conjugacy class by itself
                    self.is_inversion_symmetry_determinable = True
                    self.example_properness_conjugacy_class_pair.append(cc[0])
                sanitized_conjugacy_classes[cc[0]] = cc
                self.conjugacy_class_time_reversal[cc[0]] = rep_spatial[2]
                continue
            # First, check if cc is proper (this is uniquely determined by any element inside)
            if rep_spatial[3]:
                # improper
                # 2-fold rotoinversions are mirrors - we denote such classes as "m"
                # larger foldednesses are denoted by S_f
                if rep_spatial[1][1] == 2:
                    # mirror
                    # first, we describe the axis
                    axis_desc = "m"
                    if rep_spatial[0][0] == 1.0:
                        axis_desc += "_x"
                    elif rep_spatial[0][1] == 1.0:
                        axis_desc += "_y"
                    elif rep_spatial[0][2] == 1.0:
                        axis_desc += "_z"
                    else:
                        axis_desc += get_nondescript_desc("m")
                else:
                    # rotoinversion
                    axis_desc = "S"
                    if rep_spatial[0][0] == 1.0:
                        axis_desc += "x"
                    elif rep_spatial[0][1] == 1.0:
                        axis_desc += "y"
                    elif rep_spatial[0][2] == 1.0:
                        axis_desc += "z"
                    else:
                        axis_desc += get_nondescript_desc(f"S_{rep_spatial[1][1]}^{rep_spatial[1][0]}")
                    if rep_spatial[1][0] > 1:
                        axis_desc += f"_{rep_spatial[1][1]}^{rep_spatial[1][0]}"
                    else:
                        axis_desc += f"_{rep_spatial[1][1]}"
            else:
                # proper
                axis_desc = "C"
                if rep_spatial[0][0] == 1.0:
                    axis_desc += "x"
                elif rep_spatial[0][1] == 1.0:
                    axis_desc += "y"
                elif rep_spatial[0][2] == 1.0:
                    axis_desc += "z"
                else:
                    axis_desc += get_nondescript_desc(f"C_{rep_spatial[1][1]}^{rep_spatial[1][0]}")
                if rep_spatial[1][0] > 1:
                    axis_desc += f"_{rep_spatial[1][1]}^{rep_spatial[1][0]}"
                else:
                    axis_desc += f"_{rep_spatial[1][1]}"
                        
                    
            if self.is_double:
                if is_contracted(cc):
                    if size > 3:
                        final_name = f"{int(size / 2)}{axis_desc}+{int(size / 2)}R{axis_desc}"
                    else:
                        final_name = f"{axis_desc}+R{axis_desc}"
                    sanitized_conjugacy_classes[final_name] = cc
                    self.conjugacy_class_time_reversal[final_name] = False # this doesn't matter - all half-integer j reps will be zero here. This is just convention.
                    continue
                else:
                    reversed_cc_index = cc_i + 1 + find_reversed_class(cc, conjugacy_class_proposal[cc_i + 1:])
                    # One of the pair is time-reversed; which one it is is determined by the characters of the spinor irrep
                    skip_classes.append(reversed_cc_index)
                    if size > 1:
                        size_prefix = f"{size}"
                    else:
                        size_prefix = ""
                    time_reversal_pairs[f"{size_prefix}{axis_desc}"] = f"{size_prefix}R{axis_desc}"
                    sanitized_conjugacy_classes[f"{size_prefix}{axis_desc}"] = cc
                    sanitized_conjugacy_classes[f"{size_prefix}R{axis_desc}"] = conjugacy_class_proposal[reversed_cc_index]
                    self.conjugacy_class_time_reversal[f"{size_prefix}{axis_desc}"] = False
                    self.conjugacy_class_time_reversal[f"{size_prefix}R{axis_desc}"] = True
                    continue
            else:
                if size > 1:
                    final_name = f"{size}{axis_desc}"
                else:
                    final_name = f"{axis_desc}"
                sanitized_conjugacy_classes[final_name] = cc
                self.conjugacy_class_time_reversal[final_name] = False
                continue
        
        
        if self.is_double:
            # we'll need the spinor rep to determine which ccs are time-reversed
            
            # the spinor rep E_1/2_g (we choose the gerade version because inversion commutes w/ everything):
            #   char(E) = 2
            #   char(Rg) = -char(g)
            #   char(g) = real
            #   char(g) is nonzero for rotations with multiplicity different to [1, 2]
            #   this determines it uniquely
            
            irreps = spgrep.irreps.enumerate_unitary_irreps_from_regular_representation(self.regular_representation)
            for i in range(len(irreps)):
                irreps[i] = np.round(irreps[i], decimals = 12)
            ccs = list(sanitized_conjugacy_classes.keys())
            rep_el = []
            for i in range(len(ccs)):
                rep_el.append(self.group_elements.index(sanitized_conjugacy_classes[ccs[i]][0]))
            
            for irrep in irreps:
                is_spin_rep = True
                for i in range(len(ccs)):
                    if ccs[i] == "E":
                        if np.trace(irrep[rep_el[i]]) != 2.0:
                            is_spin_rep = False
                            break
                    if self.is_inversion_symmetry_determinable:
                        if ccs[i] == self.example_properness_conjugacy_class_pair[1]:
                            if np.trace(irrep[rep_el[i]]) != 2.0:
                                is_spin_rep = False
                                break
                    if np.trace(irrep[rep_el[i]]) == 0.0:
                        if self.element_spatial_properties[self.group_elements[rep_el[i]]][1] != [1, 2]:
                            is_spin_rep = False
                            break
                    if ccs[i] in time_reversal_pairs.keys():
                        if np.trace(irrep[rep_el[i]]) != -np.trace(irrep[rep_el[ccs.index(time_reversal_pairs[ccs[i]])]]):
                            is_spin_rep = False
                            break
                if is_spin_rep:
                    spin_rep = irrep
                    break
            if not is_spin_rep:
                print("ERROR!!! Spin rep (gerade) not found!")
                return(-1)
            
            # TODO the identifying of E_1/2 has an issue - E_1/2 must have a positive character for each non-time reversed cc. But this contains circular logic. See D3h
            # maybe this is fine, actually - if we exchange both the conjugacy classes and the irreps (E_1/2 and E_5/2), the group behaves in the same way?
            
            """print("----------SPINOR REP SELECTED-------------")
            for i in range(len(ccs)):
                cur_character = np.trace(spin_rep[rep_el[i]])
                print(ccs[i], ":", cur_character)"""
            
            self.conjugacy_class_associated_angle = {}
            for i in range(len(ccs)):
                cur_character = np.trace(spin_rep[rep_el[i]])
                cur_angle = 2.0 * np.arccos(cur_character / 2.0) # this is in [0, 2pi]
                if cur_angle >= np.pi: # time reversal: we can contract this
                    cur_angle = 2.0 * np.pi - cur_angle # this is so that char(angle) = -char(R . angle)
                self.conjugacy_class_associated_angle[ccs[i]] = cur_angle
                if ccs[i] in time_reversal_pairs.keys():
                    if cur_character < 0:
                        # we need to switch the content of ccs
                        switch_help = sanitized_conjugacy_classes[ccs[i]]
                        sanitized_conjugacy_classes[ccs[i]] = sanitized_conjugacy_classes[time_reversal_pairs[ccs[i]]]
                        sanitized_conjugacy_classes[time_reversal_pairs[ccs[i]]] = switch_help
        
        self.set_conjugacy_classes(sanitized_conjugacy_classes)
        if self.is_double:
            return(irreps) # optimalization
        else:
            return(0) # job well done
                
    
    
    def set_irreducible_representations(self, arg1, arg2 = -1):
        
        # Full specification: dict {"name" : ndarray[matrix 1, matrix 2, ...]}
        # Minimal specification: ([dimensions], [names]). if no names, names becomes labelled minimally according to altmann
        
        if type(arg1) == dict:
            self.irreps = arg1.copy()
            self.irrep_names = list(self.irreps.keys())
            self.irrep_dimensions = {}
            for irrep in self.irrep_names:
                self.irrep_dimensions[irrep] = self.irreps[irrep].shape[1]
        
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
        if self.irrep_dimensions == {}:
            if type(character_table) == dict:
                self.irrep_names = list(character_table.keys())
                for irrep in self.irrep_names:
                    self.irrep_dimensions[irrep] = int(character_table[irrep][self.conjugacy_class_names[0]])
            else:
                for i in range(len(character_table)):
                    self.irrep_dimensions[self.irrep_names[i]] = int(character_table[i][0])
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
        
        complex_conjugate_irrep_indices = {} # {index 1 : index 2}
        
        characters = []
        for irrep in irreps:
            characters.append(spgrep.representation.get_character(irrep))
        
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
        
        def classify_irreps_by_dim(set_of_irrep_indices):
            irreps_by_dim = {} # {dim1 : [[irrep1 index], [irrep2 i, irrep2 c.c. i]...], dim2 : [[irrep3], [irrep4]...]], ...}
            for i in set_of_irrep_indices:
                cur_dim = irreps[i].shape[1]
                if cur_dim in irreps_by_dim.keys():
                    # check if there's a complex conjugate irrep present
                    cc_found = False
                    for k in range(len(irreps_by_dim[cur_dim])):
                        if len(irreps_by_dim[cur_dim][k]) > 1:
                            continue
                        if check_if_complex_conjugate(i, irreps_by_dim[cur_dim][k][0]):
                            cc_found = True
                            irreps_by_dim[cur_dim][k].append(i)
                            break
                    if not cc_found:
                        irreps_by_dim[cur_dim].append([i])
                else:
                    irreps_by_dim[cur_dim] = [[i]]
            return(irreps_by_dim)
        
        def name_nonspin_irreps_classified_by_dim(irreps_by_dim, improperness_suffix = ""):
            # if improper group, specify suffix as either "_g" or "_u"
            
            final_names = []
            reorder_indices = []
        
            count_number = [[0, 0], 0, 0, 0, 0, 0]
            
            for cur_dim in sorted(irreps_by_dim.keys()):
                for j in range(len(irreps_by_dim[cur_dim])):
                    if len(irreps_by_dim[cur_dim][j]) == 1:
                        cur_i = irreps_by_dim[cur_dim][j][0]
                        reorder_indices.append(cur_i)
                        if cur_dim == 1:
                            # check for antisymmetricity ALTMANN P 63
                            if self.z_rot_elements == []:
                                count_number[cur_dim - 1][0] += 1
                                count_number[cur_dim - 1][1] += 1
                                cur_letter = f"{Group.irrep_letter_dictionary[convention][cur_dim-1][0]}/{Group.irrep_letter_dictionary[convention][cur_dim-1][1]}"
                                final_names.append(cur_letter + "_" + str(count_number[cur_dim - 1][0]) + improperness_suffix)
                            else:
                                if check_if_symmetric(cur_i):
                                    count_number[cur_dim - 1][0] += 1
                                    final_names.append(Group.irrep_letter_dictionary[convention][cur_dim-1][0] + "_" + str(count_number[cur_dim - 1][0]) + improperness_suffix)
                                else:
                                    count_number[cur_dim - 1][1] += 1
                                    final_names.append(Group.irrep_letter_dictionary[convention][cur_dim-1][1] + "_" + str(count_number[cur_dim - 1][1]) + improperness_suffix)
                        else:
                            count_number[cur_dim - 1] += 1
                            final_names.append(Group.irrep_letter_dictionary[convention][cur_dim-1] + "_" + str(count_number[cur_dim - 1]) + improperness_suffix)
                    elif len(irreps_by_dim[cur_dim][j]) == 2:
                        cur_i1 = irreps_by_dim[cur_dim][j][0]
                        cur_i2 = irreps_by_dim[cur_dim][j][1]
                        reorder_indices.append(cur_i1)
                        reorder_indices.append(cur_i2)
                        final_names.append("1" + Group.irrep_letter_dictionary[convention][cur_dim*2-1] + improperness_suffix)
                        final_names.append("2" + Group.irrep_letter_dictionary[convention][cur_dim*2-1] + improperness_suffix)
                        complex_conjugate_irrep_indices[cur_i1] = cur_i2
            return(final_names, reorder_indices)
        
        def decompose_angular_rep(j, symmetry = "g"):
            present_irreps = []
            coefs = np.zeros(len(irreps), dtype=np.complex_)
            rep = np.zeros(len(self.group_elements), dtype=np.complex_)
            for i in range(len(self.group_elements)):
                cur_angle = self.conjugacy_class_associated_angle[self.element_conjugacy_classes[self.group_elements[i]]]
                if cur_angle == 0:
                    # identity class or time reversal class - this gives (j + 1/2) / (1/2) = 2j + 1
                    rep[i] = 2 * j + 1
                else:
                    rep[i] = np.round(np.sin((j + 0.5) * cur_angle) / np.sin(0.5 * cur_angle), decimals = 8)
                # time reversal classes - for half-integer j, we flip the signs
                if self.is_double and (int(2.0 * j) % 2 == 1):
                    if self.conjugacy_class_time_reversal[self.element_conjugacy_classes[self.group_elements[i]]]:
                    #if self.element_spatial_properties[self.group_elements[i]][2] and (int(2.0 * j) % 2 == 1):
                        rep[i] *= -1
                if self.is_inversion_symmetry_determinable:
                    if self.element_spatial_properties[self.group_elements[i]][3] and symmetry == "u":
                        rep[i] *= -1
            
            """DEBUG_rep = {}
            for cc in self.conjugacy_class_names:
                DEBUG_rep[cc] = rep[self.indices_of_representative_elements[cc]]
            
            print(DEBUG_rep)"""
            
            for i in range(len(irreps)):
                for j in range(len(self.group_elements)):
                    coefs[i] += rep[j] * characters[i][j]
                if np.round(np.real(coefs[i]), decimals = 5) != 0.0:
                    present_irreps.append(i)
            return(present_irreps)
            
                
        
        
        def name_spin_irreps_classified_by_dim(irreps_by_dim, improperness_suffix = ""):
            final_names, reorder_indices = name_nonspin_irreps_classified_by_dim(irreps_by_dim, improperness_suffix)
            unlabelled_irreps = reorder_indices.copy()
            
            #print(final_names, reorder_indices)
            
            if improperness_suffix == "_u":
                symmetry = "u"
            else:
                symmetry = "g"
            for j_double in range(len(unlabelled_irreps)):
                if len(unlabelled_irreps) == 0:
                    break
                new_present_irreps = decompose_angular_rep((2 * j_double + 1) / 2.0, symmetry)
                #print(f"j = {2*j_double + 1}/2: {new_present_irreps}")
                for i in new_present_irreps:
                    if i in unlabelled_irreps:
                        name_index = reorder_indices.index(i)
                        unlabelled_irreps.remove(i)
                        final_names[name_index] += f"(j={int(2 * j_double + 1)}/2)"
            
            return(final_names, reorder_indices)
            
        
        # First classify irreps by whether or not they're spin irreps and by their gerade/ungerade properties
        if self.is_inversion_symmetry_determinable:
            noninversion_index = self.indices_of_representative_elements[self.example_properness_conjugacy_class_pair[0]]
            inversion_index = self.indices_of_representative_elements[self.example_properness_conjugacy_class_pair[1]]
        if self.is_double:
            E_index = self.indices_of_representative_elements["E"]
            R_index = self.indices_of_representative_elements["R"]
            if self.is_inversion_symmetry_determinable:
                nonspin_irreps_g = []
                nonspin_irreps_u = []
                spin_irreps_g = []
                spin_irreps_u = []
                for i in range(len(irreps)):
                    if characters[i][E_index] == -characters[i][R_index]:
                        if characters[i][noninversion_index] == characters[i][inversion_index]:
                            spin_irreps_g.append(i)
                        else:
                            spin_irreps_u.append(i)
                    else:
                        if characters[i][noninversion_index] == characters[i][inversion_index]:
                            nonspin_irreps_g.append(i)
                        else:
                            nonspin_irreps_u.append(i)
                # Here we name and classify indices for an improper double group
                nonspin_irreps_g_by_dim = classify_irreps_by_dim(nonspin_irreps_g)
                nonspin_irreps_u_by_dim = classify_irreps_by_dim(nonspin_irreps_u)
                spin_irreps_g_by_dim = classify_irreps_by_dim(spin_irreps_g)
                spin_irreps_u_by_dim = classify_irreps_by_dim(spin_irreps_u)
                
                nonspin_names_g, nonspin_reorderings_g = name_nonspin_irreps_classified_by_dim(nonspin_irreps_g_by_dim, "_g")
                nonspin_names_u, nonspin_reorderings_u = name_nonspin_irreps_classified_by_dim(nonspin_irreps_u_by_dim, "_u")
                spin_names_g, spin_reorderings_g = name_spin_irreps_classified_by_dim(spin_irreps_g_by_dim, "_g")
                spin_names_u, spin_reorderings_u = name_spin_irreps_classified_by_dim(spin_irreps_u_by_dim, "_u")
                names = nonspin_names_g + nonspin_names_u + spin_names_g + spin_names_u
                reorderings = nonspin_reorderings_g + nonspin_reorderings_u + spin_reorderings_g + spin_reorderings_u
            else:
                nonspin_irreps = []
                spin_irreps = []
                for i in range(len(irreps)):
                    if characters[i][E_index] == -characters[i][R_index]:
                        spin_irreps.append(i)
                    else:
                        nonspin_irreps.append(i)
                # Here we name and classify indices for a proper double group
                
                
                nonspin_irreps_by_dim = classify_irreps_by_dim(nonspin_irreps)
                spin_irreps_by_dim = classify_irreps_by_dim(spin_irreps)
                
                nonspin_names, nonspin_reorderings = name_nonspin_irreps_classified_by_dim(nonspin_irreps_by_dim, "")
                spin_names, spin_reorderings = name_spin_irreps_classified_by_dim(spin_irreps_by_dim, "")
                names = nonspin_names + spin_names
                reorderings = nonspin_reorderings + spin_reorderings
                
        else:
            if self.is_inversion_symmetry_determinable:
                irreps_g = []
                irreps_u = []
                for i in range(len(irreps)):
                    if characters[i][noninversion_index] == characters[i][inversion_index]:
                        irreps_g.append(i)
                    else:
                        irreps_u.append(i)
                # Here we name and classify indices for an improper group
                irreps_g_by_dim = classify_irreps_by_dim(irreps_g)
                irreps_u_by_dim = classify_irreps_by_dim(irreps_u)
                
                names_g, reorderings_g = name_nonspin_irreps_classified_by_dim(irreps_g_by_dim, "_g")
                names_u, reorderings_u = name_nonspin_irreps_classified_by_dim(irreps_u_by_dim, "_u")
                names = names_g + names_u
                reorderings = reorderings_g + reorderings_u
            else:
                # Here we name and classify indices for a proper group
                irreps_by_dim = classify_irreps_by_dim(irreps)
                names, reorderings = name_nonspin_irreps_classified_by_dim(irreps_by_dim, "")
        
        # reorder irreps, since dict remembers insertion order
        irrep_dict = {}
        for i in range(len(reorderings)):
            j = reorderings[i]
            irrep_dict[names[i]] = irreps[j]
            #irrep_dict[names[j]] = irreps[j]
        
        # find identity irrep
        self.identity_irrep
        for irrep in irrep_dict.keys():
            is_identity_irrep = True
            for i in range(len(irrep_dict[irrep])):
                if np.round(np.trace(irrep_dict[irrep][i]), decimals = 3) != 1:
                    is_identity_irrep = False
                    break
            if is_identity_irrep:
                self.identity_irrep = irrep
                break
        
        # Find the inversion irrep
        for irrep in irrep_dict.keys():
            is_inversion_irrep = True
            for i in range(len(irrep_dict[irrep])):
                if self.element_spatial_properties[self.group_elements[i]][3]:
                    # rotoinversion
                    if np.round(np.trace(irrep_dict[irrep][i]), decimals = 3) != -1:
                        is_inversion_irrep = False
                        break
                else:
                    #rotation
                    if np.round(np.trace(irrep_dict[irrep][i]), decimals = 3) != 1:
                        is_inversion_irrep = False
                        break
            if is_inversion_irrep:
                self.inversion_irrep = irrep
                break
        
        # Commit complex conjugate pairs to memory
        for i1, i2 in complex_conjugate_irrep_indices.items():
            #self.complex_conjugate_irreps[names[reorderings[i1]]] = names[reorderings[i2]]
            self.complex_conjugate_irreps[names[reorderings.index(i1)]] = names[reorderings.index(i2)]
        
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
        
        # now, we check which elements are rotations around the z axis, and commit spatial properties to memory:
        self.z_rot_elements = []
        self.element_spatial_properties = {}
        for group_element in self.group_elements:
            #if np.all(np.matmul(self.group_operations[group_element].cartesian_rep(), np.array([0.0, 0.0, 1.0])) == np.array([0.0, 0.0, 1.0])):
            if np.all(self.group_operations[group_element].axis == [0.0, 0.0, 1.0]) and self.group_operations[group_element].inversion == False:
                self.z_rot_elements.append(group_element)
            if self.group_operations[group_element].inversion:
                self.is_proper = False
            self.element_spatial_properties[group_element] = [self.group_operations[group_element].axis, self.group_operations[group_element].multiplicity, False, self.group_operations[group_element].inversion]
        
        self.order = len(self.group_elements)
        

    def get_double_group(self):
        
        # REQUIREMENTS: self.group_elements, self.group_operations, self.z_rot_elements
        # SPECIFIED PARAMETERS: self.group_elements, self.group_operations, self.z_rot_elements, self.multiplication_table
        
        # This is the group doubling method. A few notes:
        #   -This method overwrites self.group_elements, self.group_operations, and self.multiplication_table, so dependent methods have to be re-called
        #   -The resulting double group's faithful rep is no longer expressable as improper rotations, but rather as SU(2) matrices + inversions.
        
        element_matrices = [0] * self.order * 2
        
        for i in range(self.order):
            element_matrices[i]     = self.group_operations[self.group_elements[i]].SU2_rep()#ImproperSpinRotation_rep()
            element_matrices[i + self.order] = element_matrices[i].reverse()#- self.group_operations[self.group_elements[i]].SU2_rep()
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
                product_m = element_matrices[i] * element_matrices[j]#np.matmul(element_matrices[i], element_matrices[j])
                product_index = product_m.find_closest_operation(element_matrices)
                new_mt[i].append(self.group_elements[product_index])
        
        self.multiplication_table = new_mt
        # update the group operations and spatial properties
        # the multiplicity of a time reversal rotation is unconstrained - angle is bigger than 2pi
        for i in range(len(self.group_elements)):
            self.group_operations[self.group_elements[i]] = element_matrices[i]
        for i in range(self.order):
            self.element_spatial_properties[self.group_elements[i + self.order]] = self.element_spatial_properties[self.group_elements[i]].copy()
            self.element_spatial_properties[self.group_elements[i + self.order]][2] = True # we change the time reversal property
        
        # Now we update the group order and self descriptors
        self.order = len(self.group_elements)
        self.is_double = True
            

    
    
    def conjugacy_classes_from_multiplication_table(self):
        # REQUIREMENTS: self.group_elements, self.multiplication_table
        # SPECIFIED PARAMETERS: self.conjugacy_classes, self.conjugacy_class_names, self.conjugacy_class_sizes
        # this function can be called INSTEAD of self.set_conjugacy_classes()
        
        self.order = len(self.group_elements)
        # elements i, j are conjugate if there exists k such that multi[i][k] = multi[k][j]
        
        # For optimalization, we classify the elements by their order:
        element_order = {}
        orders_list = []
        for i in range(self.order): # since the maximum allowed order is the group order, which needs to be the highest ALLOWED index - so we shift indices by 1
            orders_list.append([])
        
        
        for i in range(self.order):
            j = 1
            x = i
            while(x != 0):
                j += 1
                x = self.group_elements.index(self.multiplication_table[x][i])
                
                #safety break
                if j > self.order:
                    print("ERROR: something broke in the order calculation, chief.")
                    return(-1)
            element_order[self.group_elements[i]] = j
            orders_list[j - 1].append(i)
        
        # Within each order list, categorize by classes
        
        def is_in_class_with(x, y):
            for k in range(self.order):
                if self.multiplication_table[x][k] == self.multiplication_table[k][y]:
                    return(True)
            return(False)
        
        conjugacy_classes_dict = {}
        
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
        
        # We pass down the irreps if needed
        return(self.classify_conjugacy_classes(conjugacy_classes_dict))
    
    def get_character_table(self, tolerance_decimals = 5):
        
        # REQUIREMENTS: self.irreps, self.irrep_names, self.group_elements, self.conjugacy_class_names, self.indices_of_representative_elements
        # SPECIFIED PARAMETERS: self.character_table
        # this function can be called INSTEAD of self.set_character_table()
        # to have the indices of representative elements, this function should only be called after self.conjugacy_classes_from_multiplication_table
        
        #print(self.irrep_names)
        
        characters = []
        for irrep in self.irrep_names:
            characters.append(spgrep.representation.get_character(self.irreps[irrep]))
        
        #conjugacy_class_sizes = []
        #for i in range(len(conjugacy_class_names)):
        #    conjugacy_class_sizes.append(len(conjugacy_classes_dict[conjugacy_class_names[i]]))
        
        self.character_table = np.zeros((len(self.irreps), len(self.irreps)), dtype=np.complex_)
        self.irrep_characters = {}
        
        for i_irrep in range(len(self.irreps)):
            for i_cc in range(len(self.conjugacy_class_names)):
                new_char = np.round(characters[i_irrep][self.indices_of_representative_elements[self.conjugacy_class_names[i_cc]]], decimals = tolerance_decimals)
                self.character_table[i_irrep][i_cc] = new_char
            self.irrep_characters[self.irrep_names[i_irrep]] = Representation(self, self.character_table[i_irrep])
    
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
    
    
    def find_cartesian_basis_irreps(self):
        
        # This method finds out which irreps have (x,y,z) or subsets of thereof as their basis functions
        
        # First we find out how the three cartesian unit vectors transform
        # When the 3d cartesian rep gets multiplied from the right by the identity matrix, we get a matrix of column vectors
        
        def old_is_this_basis_to_irrep(irrep, basis):
            # basis is an ORTHONORMAL set of cartesian vectors
            # D(R)_ab = <a|R|b> (Dresselhaus, p. 77)
            irrep_dim = self.irrep_dimensions[irrep]
            for i in range(len(self.group_elements)):
                # we determine the transformed basis vectors
                cur_spatial_props = self.element_spatial_properties[self.group_elements[i]]
                cur_matrix = ImproperRotation(cur_spatial_props[0], cur_spatial_props[1], cur_spatial_props[3]).cartesian_rep()
                transformed_basis = []
                for vec in basis:
                    transformed_basis.append(np.reshape(np.array(np.matmul(cur_matrix, vec)), len(vec)))
                #print(self.irreps)
                for a in range(irrep_dim):
                    for b in range(irrep_dim):
                        matrix_element = np.dot(basis[a], transformed_basis[b])
                        
                        if irrep == "E_1" and basis[0][0] == 1.0 and basis[1][1] == 1.0:
                            print(f"a:{a},b:{b} : {matrix_element} vs {self.irreps[irrep][i][a][b]}")
                        
                        
                        # NOTE this WOULD work for C3v if we took the real part - but why do we have to do that???? Why are non-spin irreps complex in double groups??
                        if np.round(matrix_element, decimals = 5) != np.round(self.irreps[irrep][i][a][b], decimals = 5):
                            return(False)
            return(True)
        
        def is_this_basis_to_irrep(irrep, basis):
            # basis is an ORTHONORMAL set of cartesian vectors
            # D(R)_ab = <a|R|b> (Dresselhaus, p. 77)
            irrep_dim = self.irrep_dimensions[irrep]
            for i in range(len(self.group_elements)):
                # we determine the transformed basis vectors
                cur_spatial_props = self.element_spatial_properties[self.group_elements[i]]
                cur_matrix = ImproperRotation(cur_spatial_props[0], cur_spatial_props[1], cur_spatial_props[3]).cartesian_rep()
                transformed_basis = []
                for vec in basis:
                    transformed_basis.append(np.reshape(np.array(np.matmul(cur_matrix, vec)), len(vec)))
                # Checks only diagonally
                current_trace = 0.0
                for a in range(irrep_dim):
                    current_trace += np.dot(basis[a], transformed_basis[a])
                # NOTE this WOULD work for C3v if we took the real part - but why do we have to do that???? Why are non-spin irreps complex in double groups??
                if np.round(current_trace, decimals = 5) != np.round(np.trace(self.irreps[irrep][i]), decimals = 5):
                    return(False)
            return(True)
            
        
        x_transformed = {}
        y_transformed = {}
        z_transformed = {}
        
        cartesian_rep = {}
        
        for cc in self.conjugacy_class_names:
            cur_spatial_props = self.element_spatial_properties[self.group_elements[self.indices_of_representative_elements[cc]]]
            cur_matrix = np.transpose(ImproperRotation(cur_spatial_props[0], cur_spatial_props[1], cur_spatial_props[3]).cartesian_rep()) #(x,y,z) cannot be basis funcs of spinor reps
            cartesian_rep[cc] = np.trace(cur_matrix)
        
        cartesian_rep_decomposition = self.reduce_representation(cartesian_rep)
        
        cartesian_unit_vectors = {"x" : np.array([1.0, 0.0, 0.0]), "y" : np.array([0.0, 1.0, 0.0]), "z" : np.array([0.0, 0.0, 1.0])}
        possible_basis_combinations = [
            [["x"]    , ["y"]    , ["z"]    ], # 1-dimensional
            [["x","y"], ["y","z"], ["z","x"]],
            [["x","y","z"]]
        ]
        self.cartesian_basis = {} # {"irrep" : ["polarisation label"]}
        for i in range(len(cartesian_rep_decomposition[0])):
            if int(cartesian_rep_decomposition[0][i]) > 0:
                self.cartesian_basis[self.irrep_names[i]] = []
                # this irrep has at least one of the three unit vectors as its basis
                cur_irrep_dim = self.irrep_dimensions[self.irrep_names[i]]
                for possible_basis in possible_basis_combinations[cur_irrep_dim-1]:
                    cur_basis_vectors = []
                    for label in possible_basis:
                        cur_basis_vectors.append(cartesian_unit_vectors[label])
                    if is_this_basis_to_irrep(self.irrep_names[i], cur_basis_vectors):
                        self.cartesian_basis[self.irrep_names[i]].append(",".join(possible_basis))
                
        """
        for i in range(len(self.group_elements)):
            element = self.group_elements[i]
            cur_spatial_props = self.element_spatial_properties[element]
            # a matrix is a list of rows, so for a list of columns we take the transpose:
            cur_matrix = np.transpose(ImproperRotation(cur_spatial_props[0], cur_spatial_props[1], cur_spatial_props[3]).cartesian_rep()) #(x,y,z) cannot be basis funcs of spinor reps
            x_transformed[element] = cur_matrix[0]
            y_transformed[element] = cur_matrix[1]
            z_transformed[element] = cur_matrix[2]
            """
            # for each irrep, we check if it could take any subset of (x,y,z) as its bases.
            # For 1D irreps, we check individual x, y, z, for 2D irreps we check (x,y), (x, z) and (y, z)
    
    
    def find_wigner_representation(self, j):
        result = {}
        for g in self.group_elements:
            result[g] = self.group_operations[g].wigner_d_matrix(j)
        return(result)
    
    
    # -------------------- initializer wrappers
    
    
    def initialize_from_multiplication_table(self):
        
        # REQUIREMENTS: self.group_elements, self.multiplication_table
        # SPECIFIED PARAMETERS: everything
        
        self.order = len(self.group_elements)
        
        # Here we find the regular representation
        self.get_regular_representation()
        # First, we find the conjugacy classes
        irreps = self.conjugacy_classes_from_multiplication_table()
        
        #for i in range(len(group_elements)):
        #    print(group_elements[i], regular_rep[i])
        
        # If irreps weren't found in conjugacy class naming, we create them now
        if not self.is_double:
            irreps = spgrep.irreps.enumerate_unitary_irreps_from_regular_representation(self.regular_representation)
            for i in range(len(irreps)):
                irreps[i] = np.round(irreps[i], decimals = 12)
        
        # We name and classify irreps
        self.name_and_set_irreps(irreps)
                
        #irrep_reordering, irrep_names = name_irreps(irreps, group_elements, rotation_matrices)
        #irreps = [irreps[i] for i in irrep_reordering]
        
        # We set the character table
        self.get_character_table()
        
        # Determine cartesian basis
        self.find_cartesian_basis_irreps()
        
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
        
        
    
    
    # ------------------------ representation methods
    
    def reduce_representation(self, reducible_representation):
        
        # This is either a dict {"conjugacy class name" : character} or a list [cc_i] = character
        
        coefs = [0] * len(self.character_table)
        
        for i in range(len(self.character_table)):
            for j in range(len(self.character_table[i])):
                if type(reducible_representation) == Representation:
                    coefs[i] += reducible_representation.characters[self.conjugacy_class_names[j]] * self.character_table[i][j] * self.conjugacy_class_sizes[j]
                elif type(reducible_representation) == dict:
                    #if self.conjugacy_class_names[j] in reducible_representation:
                    coefs[i] += reducible_representation[self.conjugacy_class_names[j]] * self.character_table[i][j] * self.conjugacy_class_sizes[j]
                else:
                    coefs[i] += reducible_representation[j] * self.character_table[i][j] * self.conjugacy_class_sizes[j]
            coefs[i] /= self.order
            coefs[i] = np.round(coefs[i], decimals = 3)
            
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
    
    def angular_representation(self, j, inversion_symmetry = "unspecified"):
        # for a given j value, this creates the reducible angular representation (as a subgroup of the full rotation group)
        # if symmetry = "g" (gerade), then characters dont flip sign for inversions; for symmetry = "u" (ungerade), they do
        
        
        # This is the character of an irrep based on the Wigner D-matrices.
        # For odd j = ungerade, for even j = gerade. For half-integer j: no idea but let's use gerade as a convention
        # Idea: with (-1)^j, half-integer js should multiply rotoinversions by i or -i (1/2 or 3/2). Try this!!
        
        if inversion_symmetry == "unspecified":
            if np.round(j) == j:
                if j % 2 == 0: #evem
                    inversion_symmetry = "g"
                elif j % 2 == 1: #odd
                    inversion_symmetry = "u"
        
        # One way to resolve the half-integer j case is to say: from C_i we see that inversion commutes with R, so we can say
        # j=1/2 has a special rep (which is the E(j=1/2) irrep) and then the angular rep of l+1/2 = ang.rep(l) + E(j=1/2). Yay!!
        
        #inversion_factor = np.power(1j, j*2)
        
        rep = {}
        is_j_half_integer = (int(2.0 * j) % 2 == 1)
        for cc in self.conjugacy_class_names:
            cur_mult = self.element_spatial_properties[self.group_elements[self.indices_of_representative_elements[cc]]][1]
            cur_angle = self.conjugacy_class_associated_angle[cc]
            if cur_angle == 0:
                # identity class or time reversal class - this gives (j + 1/2) / (1/2) = 2j + 1
                rep[cc] = 2 * j + 1
            else:
                rep[cc] = np.round(np.sin((j + 0.5) * cur_angle) / np.sin(0.5 * cur_angle), decimals = 8)
            # time reversal classes - for half-integer j, we flip the signs
            if self.is_double:
                if self.conjugacy_class_time_reversal[cc] and is_j_half_integer:
                    rep[cc] *= -1
            #if self.is_inversion_symmetry_determinable:
            if self.element_spatial_properties[self.group_elements[self.indices_of_representative_elements[cc]]][3] and inversion_symmetry == "u":
                rep[cc] *= -1
                
        return(Representation(self, rep))
    
    def does_rep_contain_identity(self, rep):
        # checks whether representation contains the identity representation (which has characters 1 for all CCs)
        id_coef = 0
        
        for i in range(len(self.conjugacy_class_names)):
            if type(rep) == Representation:
                id_coef += rep.characters[self.conjugacy_class_names[i]] * self.conjugacy_class_sizes[i]
            elif type(rep) == dict:
                id_coef += rep[self.conjugacy_class_names[i]] * self.conjugacy_class_sizes[i]
            else:
                id_coef += rep[i] * self.conjugacy_class_sizes[i]
        id_coef /= self.order
        id_coef = np.round(np.real(id_coef))
        return(id_coef != 0.0)
    
    def separate_constituent_representations(self, set_of_irreps, clump_conjugate_irreps = True):
        # if set_of_irreps is a Representation, we reduce it. Otherwise it must be a list of coefficients.
        if type(set_of_irreps) == Representation:
            set_of_irreps = self.reduce_representation(set_of_irreps)[0]
        result = {} # {"unique_label" : Rep}
        # we assume all coeffs are NON-NEGATIVE INTEGERS
        for i in range(len(set_of_irreps)):
            if set_of_irreps[i] == 0.0:
                continue
            if clump_conjugate_irreps and self.irrep_names[i] in self.complex_conjugate_irreps.keys():
                set_of_irreps[self.irrep_names.index(self.complex_conjugate_irreps[self.irrep_names[i]])] -= int(set_of_irreps[i]) # assuming this is always the same
                for j in range(int(set_of_irreps[i])):
                    result[f"{self.irrep_names[i]}+{self.complex_conjugate_irreps[self.irrep_names[i]]}[{j+1}]"] = self.irrep_characters[self.irrep_names[i]] + self.irrep_characters[self.complex_conjugate_irreps[self.irrep_names[i]]]
            else:
                for j in range(int(set_of_irreps[i])):
                    result[f"{self.irrep_names[i]}[{j+1}]"] = self.irrep_characters[self.irrep_names[i]]#self.irrep_names[i]
        return(result)
    
    """def allowed_transitions_between_reps(self, rep1, rep2, interaction_term_rep):
        # rep1 and rep2 can be reducible - we reduce them and then treat each component as a separate energy level
        # Returns {"polarisation" : [allowed transitions, dark transitions]}
        rep1_reduction, hr1 = self.reduce_representation(rep1)
        rep2_reduction, hr2 = self.reduce_representation(rep2)
        
        energy_levels1 = self.separate_constituent_representations(rep1_reduction)
        energy_levels2 = self.separate_constituent_representations(rep2_reduction)
        
        
        
        allowed_transitions = [] # ["label1 -> label2"]
        dark_transitions = []
        for E1 in energy_levels1.keys():
            for E2 in energy_levels2.keys():
                transition_rep = self.irrep_characters[energy_levels1[E1]] * self.irrep_characters[energy_levels1[E2]] * interaction_term_rep
                if self.does_rep_contain_identity(transition_rep):
                    allowed_transitions.append(f"{E1} -> {E2}")
                else:
                    dark_transitions.append(f"{E1} -> {E2}")
        return(allowed_transitions, dark_transitions)"""
    
    def reduce_representation_and_divide_basis(self, representation, clump_conjugate_irreps = True):
        
        # This method not only reduces the rep, but also assigns the indices in the original basis (D-dim vector) to each constituent irrep 
        
        # representation is either a Representation, a list of matrices ordered the same way as self.group_elements or a dictionary with said group_elements as keys
        
        reducible_rep = []
        if type(representation) == dict:
            for g in self.group_elements:
                reducible_rep.append(representation[g])
        else:
            reducible_rep = representation.copy()
        reducible_rep = np.array(reducible_rep)
        
        D = reducible_rep.shape[1]
        
        conjugate_class_characters = {}
        for cc in self.conjugacy_class_names:
            conjugate_class_characters[cc] = np.trace(reducible_rep[self.indices_of_representative_elements[cc]])
        
        constituent_reps = list(self.separate_constituent_representations(Representation(self, conjugate_class_characters), clump_conjugate_irreps).values())
        
        unclassified_basis_indices = []
        for i in range(D):
            unclassified_basis_indices.append(i)
        
        resulting_irrep_basis = [] # [[Representation, [i1, i2, i3...]], ...]
        
        for i in range(len(constituent_reps)):
            cur_dim = int(np.real(constituent_reps[i].characters[self.example_properness_conjugacy_class_pair[0]]))
            possible_basis_indices_i = index_sublists(len(unclassified_basis_indices), cur_dim)
            basis_indices_sublist_found = False
            for j in range(len(possible_basis_indices_i)):
                possible_basis_indices = []
                for index in possible_basis_indices_i[j]:
                    possible_basis_indices.append(unclassified_basis_indices[index])
                # Check if these basis indices are correct
                is_basis_indices_subset_correct = True
                for cc in self.conjugacy_class_names:
                    cur_subtrace = 0.0
                    for b_i in possible_basis_indices:
                        cur_subtrace += reducible_rep[self.indices_of_representative_elements[cc]][b_i][b_i]
                    if np.round(cur_subtrace, decimals = 5) != np.round(constituent_reps[i].characters[cc]):
                        is_basis_indices_subset_correct = False
                        break
                if is_basis_indices_subset_correct:
                    basis_indices_sublist_found = True
                    break
            if not basis_indices_sublist_found:
                # whut?
                print(f"ERROR: subtrace index sublist not found for rep {self.reduce_representation(constituent_reps[i])[1]}")
                return(-1)
            else:
                resulting_irrep_basis.append([constituent_reps[i], possible_basis_indices])
                # We delete the used basis indices
                possible_basis_indices_i[j].reverse()
                for used_basis_index in possible_basis_indices_i[j]:
                    del unclassified_basis_indices[used_basis_index]
        return(resulting_irrep_basis)
    
    
    
    def allowed_transitions_between_reps(self, rep1, rep2, clump_conjugate_irreps = True):
        # rep1 and rep2 can be reducible - we reduce them and then treat each component as a separate energy level
        # Returns {"polarisation" : [allowed transitions, dark transitions]}
        rep1_reduction, hr1 = self.reduce_representation(rep1)
        rep2_reduction, hr2 = self.reduce_representation(rep2)
        
        energy_levels1 = self.separate_constituent_representations(rep1_reduction, clump_conjugate_irreps)
        energy_levels2 = self.separate_constituent_representations(rep2_reduction, clump_conjugate_irreps)
        
        result = {}
        for interaction_irrep in self.cartesian_basis.keys():
            # All polarisations with the same dipole term rep will have the same transitions, so we combine them
            allowed_transitions = [] # ["label1 -> label2"]
            dark_transitions = []
            for E1 in energy_levels1.keys():
                for E2 in energy_levels2.keys():
                    transition_rep = energy_levels1[E1] * energy_levels2[E2] * self.irrep_characters[interaction_irrep]
                    if self.does_rep_contain_identity(transition_rep):
                        allowed_transitions.append(f"{E1} -> {E2}")
                    else:
                        dark_transitions.append(f"{E1} -> {E2}")
            result["(" + "); (".join(self.cartesian_basis[interaction_irrep]) + ")"] = [allowed_transitions.copy(), dark_transitions.copy()]
        return(result)
    
    
    
    # ---------------------- group methods
    
    def add_subgroup(self, subgroup):
        
        # initializes subgroup relations
        # NOTE "subgroup" needs to be an instance of the class Group
        
        # REQUIREMENTS: both self and subgroup need their group_operations to contain a faithful rep (not necessarily an ImproperRotation)
        
        new_element_relations = {}
        new_conjugacy_relations = {}
        
        for s_i in range(len(subgroup.group_elements)):
            related_element_found = False
            for i in range(len(self.group_elements)):
                if np.all(subgroup.group_operations[subgroup.group_elements[s_i]] == self.group_operations[self.group_elements[i]]):
                    related_element_found = True
                    break
            if not related_element_found:
                # the input subgroup is not a true subgroup of self
                print(f"ERROR: The input subgroup {subgroup.name} is not a true subgroup of {self.name}.")
                return(-1)
            new_element_relations[subgroup.group_elements[s_i]] = self.group_elements[i]
            new_conjugacy_relations[subgroup.element_conjugacy_classes[subgroup.group_elements[s_i]]] = self.element_conjugacy_classes[self.group_elements[i]]
        
        self.subgroups.append(subgroup.name)
        self.subgroup_element_relations[subgroup.name] = [subgroup, new_element_relations]
        self.subgroup_conjugacy_relations[subgroup.name] = [subgroup, new_conjugacy_relations]
        
    def rep_to_subgroup_rep(self, subgroup_name, representation):
        
        # subgroup_name is the name of an initialized subgroup
        # representation is an instance of the Representation class
        # REQUIREMENTS: subgroup relations
        new_rep = {}
        if type(representation) == Representation:
            for sub_cc in self.subgroup_conjugacy_relations[subgroup_name][1].keys():
                new_rep[sub_cc] = representation.characters[self.subgroup_conjugacy_relations[subgroup_name][1][sub_cc]]
        elif type(representation) == dict:
            for sub_cc in self.subgroup_conjugacy_relations[subgroup_name][1].keys():
                new_rep[sub_cc] = representation[self.subgroup_conjugacy_relations[subgroup_name][1][sub_cc]]
        else:
            # a list over conjugacy classes
            for sub_cc in self.subgroup_conjugacy_relations[subgroup_name][1].keys():
                cc = self.subgroup_conjugacy_relations[subgroup_name][1][sub_cc]
                i = self.conjugacy_class_names.index(cc)
                new_rep[sub_cc] = representation[i]
        return(Representation(self.subgroup_conjugacy_relations[subgroup_name][0], new_rep))
        


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


a = Cz_4.SU2_rep() * Cz_4.SU2_rep()
b = Cz_2.SU2_rep().reverse()

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
