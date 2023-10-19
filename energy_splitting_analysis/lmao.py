
import numpy as np
from groups import Group

# char table for D4h double group

D4h = 32

D4h_conjugacy_classes = ["E", "2C4", "C2+RC2", "2C'+2RC'", "2C''+2RC''", "i", "2S4", "sigma_h+R.sigma_h", "2 sigma_v+2R.sigma_v", "2 sigma_d + 2R.sigma_d", "R", "2RC4", "Ri", "2RS4"]

D4h_cc_size = [1, 2, 2, 4, 4, 1, 2, 2, 4, 4, 1, 2, 1, 2]

D4h_irreps = ["A1g", "A2g", "B1g", "B2g", "Eg", "A1u", "A2u", "B1u", "B2u", "Eu", "E1/2g", "E3/2g", "E1/2u", "E3/2u"]

#print(len(conjugacy_classes))

D4h_char_table = [[1]*14,
              [1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, 1],
              [1, -1, 1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1],
              [1, -1, 1, -1, 1, 1, -1, 1, -1, 1, 1, -1, 1, -1],
              [2, 0, -2, 0, 0, 2, 0, -2, 0, 0, 2, 0, 2, 0],
              [1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1],
              [1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1],
              [1, -1, 1, 1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1],
              [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1],
              [2, 0, -2, 0, 0, -2, 0, 2, 0, 0, 2, 0, -2, 0],
              [2, np.sqrt(2), 0, 0, 0, 2, np.sqrt(2), 0, 0, 0, -2, -np.sqrt(2), -2, -np.sqrt(2)],
              [2, -np.sqrt(2), 0, 0, 0, 2, -np.sqrt(2), 0, 0, 0, -2, np.sqrt(2), -2, np.sqrt(2)],
              [2, np.sqrt(2), 0, 0, 0, -2, -np.sqrt(2), 0, 0, 0, -2, -np.sqrt(2), 2, np.sqrt(2)],
              [2, -np.sqrt(2), 0, 0, 0, -2, np.sqrt(2), 0, 0, 0, -2, np.sqrt(2), 2, -np.sqrt(2)]
              ]


        

D4h = Group(D4h_conjugacy_classes, D4h_irreps, D4h_cc_size, D4h_char_table)

Ih_conjugacy_classes = ["E", "12C_5", "12C_5^2", "20C_3", "15C_2+15RC_2", "i", "12S_10^3", "12S_10", "20S_6", "15.sigma+15R.sigma", "R", "12RC_5", "12RC_5^2", "20RC_3", "Ri", "12RS_10^3", "12RS_10", "20RS_6"]
Ih_cc_size = [1, 12, 12, 20, 30, 1, 12, 12, 20, 30, 1, 12, 12, 20, 1, 12, 12, 20]

Ih_irreps = ["Ag", "T1g", "T2g", "Fg", "Hg", "Au", "T1u", "T2u", "Fu", "Hu", "E1/2g", "E7/2g", "E3/2g", "E5/2g", "E1/2u", "E7/2u", "E3/2u", "E5/2u"]

def c(r, e=1.0):
    return(np.cos(np.pi * e / r))

Ih_char_table = [
    [1]*18, #Ag
    [3, 2*c(5), 2*c(5, 3), 0, -1, 3, 2*c(5), 2*c(5, 3), 0, -1, 3, 2*c(5), 2*c(5, 3), 0, 3, 2*c(5), 2*c(5, 3), 0], #T1g
    [3, 2*c(5, 3), 2*c(5), 0, -1, 3, 2*c(5, 3), 2*c(5), 0, -1, 3, 2*c(5, 3), 2*c(5), 0, 3, 2*c(5, 3), 2*c(5), 0], #T2g
    [4, -1, -1, 1, 0, 4, -1, -1, 1, 0, 4, -1, -1, 1, 4, -1, -1, 1], #Fg
    [5, 0, 0, -1, 1, 5, 0, 0, -1, 1, 5, 0, 0, -1, 5, 0, 0, -1], #Hg
    [1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1], #Au
    [3, 2*c(5), 2*c(5, 3), 0, -1, -3, -2*c(5), -2*c(5, 3), 0, 1, 3, 2*c(5), 2*c(5, 3), 0, -3, -2*c(5), -2*c(5, 3), 0], #T1u
    [3, 2*c(5, 3), 2*c(5), 0, -1, -3, -2*c(5, 3), -2*c(5), 0, 1, 3, 2*c(5, 3), 2*c(5), 0, -3, -2*c(5, 3), -2*c(5), 0], #T2u
    [4, -1, -1, 1, 0, -4, 1, 1, -1, 0, 4, -1, -1, 1, -4, 1, 1, -1], #Fu
    [5, 0, 0, -1, 1, -5, 0, 0, 1, -1, 5, 0, 0, -1, -5, 0, 0, 1], #Hu
    [2, 2*c(5), 2*c(5, 2), 1, 0, 2, 2*c(5), 2*c(5, 2), 1, 0, -2, -2*c(5), -2*c(5, 2), -1, -2, -2*c(5), -2*c(5, 2), -1], #E1/2g
    [2, 2*c(5, 3), 2*c(5, 4), 1, 0, 2, 2*c(5, 3), 2*c(5, 4), 1, 0, -2, -2*c(5, 3), -2*c(5, 4), -1, -2, -2*c(5, 3), -2*c(5, 4), -1], #E7/2g
    [4, 1, -1, -1, 0, 4, 1, -1, -1, 0, -4, -1, 1, 1, -4, -1, 1, 1], #E3/2g
    [6, -1, 1, 0, 0, 6, -1, 1, 0, 0, -6, 1, -1, 0, -6, 1, -1, 0], #E5/2g
    [2, 2*c(5), 2*c(5, 2), 1, 0, -2, -2*c(5), -2*c(5, 2), -1, 0, -2, -2*c(5), -2*c(5, 2), -1, 2, 2*c(5), 2*c(5, 2), 1], #E1/2u
    [2, 2*c(5, 3), 2*c(5, 4), 1, 0, -2, -2*c(5, 3), -2*c(5, 4), -1, 0, -2, -2*c(5, 3), -2*c(5, 4), -1, 2, 2*c(5, 3), 2*c(5, 4), 1], #E7/2u
    [4, 1, -1, -1, 0, -4, -1, 1, 1, 0, -4, -1, 1, 1, 4, 1, -1, -1], #E3/2u
    [6, -1, 1, 0, 0, -6, 1, -1, 0, 0, -6, 1, -1, 0, 6, -1, 1, 0] #E5/2u
    ]


Ih_alphas = [0, 2.0 * np.pi / 5.0, 2.0 * np.pi * 2.0 / 5.0, 2.0 * np.pi / 3.0, np.pi, 0, 2.0 * np.pi / 5.0, 2.0 * np.pi * 2.0 / 5.0, 2.0 * np.pi / 3.0, np.pi, 0, 2.0 * np.pi / 5.0, 2.0 * np.pi * 2.0 / 5.0, 2.0 * np.pi / 3.0, 0, 2.0 * np.pi / 5.0, 2.0 * np.pi * 2.0 / 5.0, 2.0 * np.pi / 3.0]

Ih_time_reversal = [1]*10 + [-1] * 8

Ih = Group(Ih_conjugacy_classes, Ih_irreps, Ih_cc_size, Ih_char_table)




D4h_reducible_rep = [16, 0, 0, 0, 0, 16, 0, 0, 0, 0, -16, 0, -16, 0]

Ih_reducible_rep = [16, -1, 1, -1, 0, 16, -1, 1, -1, 0, -16, 1, -1, 1, -16, 1, -1, 1]

"""def print_char_table(cc, ir, table):
    os = "    | "
    for i in range(len(cc)):
        os += str(cc[i]) + " | "
    print(os[:-3])
    
    for i in range(len(ir)):
        os = str(ir[i]) + " | "
        for j in range(len(table[i])):
            os += str(table[i][j]) + "   "
        print(os)"""
"""
print_char_table(D4h_conjugacy_classes, D4h_irreps, D4h_char_table)



def decompose_rep(reducible_rep, char_table, cc_size, irreps, G):

    coefs = [0] * len(char_table)

    for i in range(len(char_table)):
        #s = 0
        for j in range(len(char_table[i])):
            coefs[i] += reducible_rep[j] * char_table[i][j] * cc_size[j]
        coefs[i] /= G
        
    human_readable_output = ""
    for i in range(len(coefs)):
        if coefs[i] != 0.0:
            human_readable_output += str(int(coefs[i])) + "." + irreps[i] + " + "

    return(coefs, human_readable_output[:-3])"""

D4h_coefs, D4h_log = D4h.reduce_representation(D4h_reducible_rep)
print(D4h_log)
print("--------------------- I_h symmetry analysis --------------------------")
Ih_coefs, Ih_log = Ih.reduce_representation(Ih_reducible_rep)
print(Ih_log)


def construct_Ih_j_rep(j):
    chars = []
    for i in range(len(Ih_alphas)):
        if Ih_alphas[i] == 0:
            chars.append(2 * j + 1.0)
        else:
            chars.append(np.sin((j + 1/2.0) * Ih_alphas[i]) / np.sin(Ih_alphas[i] / 2.0))
        # if half-integer, time reversal flips spin
        if (2 * j) % 2 == 1:
           chars[i] *= Ih_time_reversal[i]
        chars[i] = np.round(chars[i], decimals = 14)
    return(chars)
            



# a j=7/2 rep
Ih_sevenhalf = [8, 2*c(5, 3)-1, 2*c(5, 4)+1, 1, 0, 8, 2*c(5, 3)-1, 2*c(5, 4)+1, 1, 0, -8, -(2*c(5, 3)-1), -(2*c(5, 4)+1), -1, -8, -(2*c(5, 3)-1), -(2*c(5, 4)+1), -1]

print(Ih_sevenhalf)
print(construct_Ih_j_rep(7/2))

print("for a j=7/2 spinor, we have the following splitting")
Ih_coefs, Ih_log = Ih.reduce_representation(Ih_sevenhalf)
print(Ih_log)

# this might be because there is no pure l=3 state (no 7-fold degeneracy in og group), so a j=7/2 state cannot be constructed as a pure state
# (which would explain why for EVERY group there is a cutoff (a maximal possible pure spinor state)

# test: if H_g labels l=2, then a spatial H_g that has a s=1/2=E_1/2 spinor part should transform as (j=5/2=I_5/2 state) + (j=3/2=E_3/2 state)

s_two = Ih_char_table[4]
s_half = Ih_char_table[10]

print(s_two, s_half)

reducible_rep = []
for i in range(len(s_two)):
    reducible_rep.append(s_two[i] * s_half[i])

Ih_coefs, Ih_log = Ih.reduce_representation(reducible_rep)
print(Ih_log)


# WHat's up with the F guy? try reducing a j=3 rep
Ih_three = construct_Ih_j_rep(3)#[7, (2*c(5, 3)-1), -2*c(5, 3), 1, -1, 7, (2*c(5, 3)-1), -2*c(5, 3), 1, -1, 7, (2*c(5, 3)-1), -2*c(5, 3), 1, 7, (2*c(5, 3)-1), -2*c(5, 3), 1]

# IT'S TRUE: j=3 state is just T2 + F



for j in range(8):

    Ih_coefs, Ih_log = Ih.reduce_representation(construct_Ih_j_rep(j))
    print(f"j = {j}:", Ih_log)
    Ih_coefs, Ih_log = Ih.reduce_representation(construct_Ih_j_rep(j+1/2))
    print(f"j = {(2*j+1)}/2:", Ih_log)
#Ih_coefs, Ih_log = Ih.reduce_representation(Ih_three)
#print("j = 3:", Ih_log)




# now, for the fun part: take the product of two reps with j1, j2. They should split into reps (NOT irreps) associated with j1-j2, j1-j2+1,...j1+j2

def decompose_aligned_j(j1, j2):
    
    rep1 = construct_Ih_j_rep(j1)
    rep2 = construct_Ih_j_rep(j2)
    
    product_rep = []
    for i in range(len(rep1)):
        product_rep.append(rep1[i] * rep2[i])
    
    #now construct the j reps
    degeneracy = int(j1 + j2 - np.abs(j1-j2) + 1)
    full_split_rep = [0]*len(product_rep)
    for i in range(degeneracy):
        j = np.abs(j1-j2) + i
        j_rep = construct_Ih_j_rep(j)
        #print(j_rep)
        for i2 in range(len(full_split_rep)):
            full_split_rep[i2] += j_rep[i2]
    
    for i in range(len(full_split_rep)):
        full_split_rep[i] = np.round(full_split_rep[i], decimals = 10)
    
    print("Gamma(j=|j1-j2|) + Gamma(j=|j1-j2|+1) + ... + Gamma(j=j1+j2):", full_split_rep)
    print("--------------------------------------------")
    print("Gamma(j1) x Gamma(j2):", product_rep)
    
    difference_magnitude = 0
    for i in range(len(product_rep)):
        difference_magnitude += (product_rep[i] - full_split_rep[i]) * (product_rep[i] - full_split_rep[i])
    difference_magnitude = np.round(difference_magnitude, decimals = 10)
    if difference_magnitude == 0:
        print("These two representations are identical!")
    else:
        print("These two representations differ by", difference_magnitude)

decompose_aligned_j(7, 7/2)








