# ---------------------------------------------------------
#
#   The leading src file of the AReTDoG package
#   AReTDoG: Algorithmic Representation Theory of Double Groups
#   Import this file to access the full scope of AReTDoG
#
#   This file provides:
#       class QDGroup(Group)
# ---------------------------------------------------------
#   Created: 8 November 2023
#   Author: Michal Horansky (michal.horansky20@imperial.ac.uk)
#
# ---------------------------------------------------------


from .class_group import *


class QDGroup(Group):
    # Instances of this class possess properties specific to the semiconductor structure, like the conduction and the valence band orbitals, the notion of holes etc
    
    def __init__(self, name):
        self.holes = {} #{"hole name" : Rep}
        self.hole_irreps = {} #{"hole name" : "hole irrep or duo of conjugate irreps"}
        self.electrons = {} # we can have multiple characters of electrons when in higher orbital (very excited)
        self.electron_irreps = {} # "rep name"
        super().__init__(name)
    
    
    
    # --------------- exciton labelling methods
    
    
    def classify_fermions(self, orbital, spin=1/2, band="valence"):
        # for GaAs: conduction band usually has s-orbital and valence band usually has p-orbital
        # if band="valence", we select highest possible j (highest energy). If "conduction", we select lowest j
        
        orbital_dict = {"s" : 0, "p" : 1, "d" : 2, "f" : 3}
        if type(orbital) == str: #otherwise must be a number
            orbital = orbital_dict[orbital]
        if band == "valence" or band == "v":
            j_val = orbital + spin
        elif band == "conduction" or band == "c":
            j_val = np.abs(orbital - spin)
        gerade_rep = self.angular_representation(j_val, "g")
        if orbital % 2 == 1:
            gerade_rep *= self.irrep_characters[self.inversion_irrep]
        
        # now we classify these:
        return(self.separate_constituent_representations(self.reduce_representation(gerade_rep)[0]))
    
    def classify_electrons(self, orbital):
        electron_dict = self.classify_fermions(orbital, 1/2, "conduction")
        electron_i = 1
        for irrep_label, rep in electron_dict.items():
            self.electrons[f"e{electron_i}"] = rep
            self.electron_irreps[f"e{electron_i}"] = self.reduce_representation(rep)[1]
            electron_i += 1
    
    def classify_holes(self, orbital):
        hole_dict = self.classify_fermions(orbital, 1/2, "valence")
        hole_i = 1
        for irrep_label, rep in hole_dict.items():
            self.holes[f"h{hole_i}"] = rep
            self.hole_irreps[f"h{hole_i}"] = self.reduce_representation(rep)[1]
            hole_i += 1
    
    
    def exciton_rep(self, electron_numbers, hole_numbers):
        result_rep = self.irrep_characters[self.irrep_names[0]] #identity rep
        i = 0
        for e_name, e_rep in self.electrons.items():
            # pauli exclusion - each filled level is labelled by identity
            if type(electron_numbers) == dict:
                free_e = electron_numbers[e_name] % int(e_rep.characters["E"])
            else:
                free_e = electron_numbers[i] % int(e_rep.characters["E"])
            
            # what to do when multiple free e? idk mate
            for j in range(free_e):
                result_rep *= e_rep
            i += 1
        i = 0
        for h_name, h_rep in self.holes.items():
            # pauli exclusion - each filled level is labelled by identity
            if type(hole_numbers) == dict:
                free_h = hole_numbers[h_name] % int(h_rep.characters["E"])
            else:
                free_h = hole_numbers[i] % int(h_rep.characters["E"])
            
            # what to do when multiple free e? idk mate
            for j in range(free_h):
                result_rep *= h_rep
            i += 1
        
        return(result_rep)
        
        
        
    
    
    
