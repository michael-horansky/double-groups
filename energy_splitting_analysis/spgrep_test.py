
import spgrep

import numpy as np

#from spgrep import get_spacegroup_irreps
#from spgrep.representation import get_character








# C3v group

print(spgrep.pointgroup.pg_dataset["3m"][0]) # this is in fractional coordinates, not cartesian coordinates!!!

rotations = np.array(spgrep.pointgroup.pg_dataset["3m"][0])
irreps = spgrep.irreps.enumerate_unitary_irreps(rotations)[0]

characters = []
for irrep in irreps:
    characters.append(spgrep.representation.get_character(irrep))
    

print(characters)

# D3h group

D3h_rots = spgrep.pointgroup.pg_dataset["-6m2"][0]

#D3h_frac_coords = [[np.cos(2.0 / 3.0 * np.pi), np.sin(2.0 / 3.0 * np.pi), 0.0], [np.cos(4.0 / 3.0 * np.pi), np.sin(4.0 / 3.0 * np.pi), 0.0], [1.0, 0.0, 0.0]]

#rotations = np.fromfunction(lambda i : np.matrix(D3h_rots[i]), (len(D3h_rots)))
rotations = np.array(D3h_rots)
irreps = spgrep.irreps.enumerate_unitary_irreps(rotations)[0]

characters = []
print("---------------------------------")
for irrep in irreps:
    characters.append(spgrep.representation.get_character(irrep))
    print(characters[-1])

new_rotations = []
for i in rotations:
    new_rotations.append(np.matrix(i))
rotations = np.array(new_rotations)

print(rotations)

sigma = rotations[3]


def find_conjugating_rotation(rots, rot1, rot2):
    
    tol = 1e-5
    
    def are_matrices_equal(a, b):
        #print(a == b)
        return(np.all(a == b))
    
    list_of_answers = []
    list_of_indices = []
    for i in range(len(rots)):
        cur_conjugator = rots[i].copy()
        cur_conj_inv = np.linalg.inv(cur_conjugator)
        if are_matrices_equal(np.matmul(np.matmul(cur_conjugator, rot1), cur_conj_inv), rot2):
            list_of_answers.append(cur_conjugator)
            list_of_indices.append(i)
    return(list_of_answers, list_of_indices)



#TODO find conjugacy classes -> then add R and see what happens
        
        

C21 = rotations[6]
C22 = rotations[7]
C23 = rotations[8]

print("C21 =", C21)
print("C22 =", C22)
print("C23 =", C23)
"""
#print(np.matmul(np.matmul(sigma, C23), sigma))
a, b = find_conjugating_rotation(rotations, C21, C22)
print("conjugators of C21 and C22:", a)"""


#print(characters)



# Constructing a double group - we have D(RA) = +- D(A), where plus and minus is for even and odd dimensional compatibility relations of the full rotation group respectively
# (so j half-integer irreps have D(A)=-D(A) and integer j irreps have D(RA)=D(A)
# However, this means that 3D spatial rotations no longer produce a faithful representation
# we need to construct the refular rep



#---------------------------------------------------------------------------------
# We embed the double group in SU(2), a double cover group of SO(3)

def SO3_to_SU2(M):
    
    # M is a list of ROWS (so that M_ij = M[i-1][j-1])
    
    # https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Rotation_matrix_%E2%86%94_Euler_axis/angle
    
    #First, we calculate the Euler angle:
    
    phi = np.arccos((M[0][0] + M[1][1] + M[2][2] - 1.0)/2.0)
    
    print(phi)
    
    return(phi)


v = np.array([np.array([0.0]), np.array([-5.0]), np.array([0.0])])

print("v turns into", np.matmul(C21, v))

print("Euler angle of C21 =", SO3_to_SU2(C21) / np.pi, "pi")



















