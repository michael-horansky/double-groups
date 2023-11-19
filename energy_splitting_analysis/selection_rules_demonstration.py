


from aretdog.class_QD_group import *




C3v = QDGroup("C3v QD")
C3v.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C3v.print_character_table()
C3v.classify_excitons(max_e = 3, max_h = 2)
C3v.print_exciton_complexes()

D3h = QDGroup("D3h QD")
D3h.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "m" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True), "m'" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], True)})
D3h.print_character_table()
D3h.classify_excitons(max_e = 3, max_h = 2)
D3h.print_exciton_complexes()

D3h.find_transition_chain()
print(D3h.transition_chain)


print("----------------- Wigner analysis -------------------")

my_j = 3/2
wigner_rep = C3v.find_wigner_representation(my_j)
print(C3v.reduce_representation(C3v.angular_representation(my_j))[1]) #NOTE add inversions to wigner d matrices - done but check if it's actually sensible pls
"""
print(wigner_rep)

for cc in C3v.conjugacy_class_names:
    print(f"{cc} : {np.trace(wigner_rep[C3v.conjugacy_classes[cc][0]])}")"""

reduced_wigner = C3v.reduce_representation_and_divide_basis(wigner_rep)
print("Wigner for j =", my_j)
for rep in reduced_wigner:
    print(f"    {C3v.reduce_representation(rep[0])[1]} : {rep[1]}")


print("Geradeness analysis")

"""
For D3h, geradeness is a problem. The j=3/2 holes are ungerade, since they're prepared from spin-orbit coupling of gamma_j=1 (the p orbit, ungerade) and gamma_j=1/2 (half-spin
particle, gerade). So, when they couple with j=1/2 electrons (gerade). the result should be wholly ungerade (in the sum). But we get a j=1 and a j=2 state - the latter would
by direct construction always be gerade! what's going on? What will it be? Is the geradeness propagating sthrough coupling in a strange manner?
"""

D3h_X10 = D3h.excitons["1X[1][1,0]"] # electron + heavy
D3h_X01 = D3h.excitons["1X[1][0,1]"] # electron + light
D3h_X = D3h_X10 + D3h_X01 # electron + hole

spherical_harmonic_prediction = D3h.angular_representation(1) + D3h.angular_representation(2)
spherical_harmonic_ungerade = D3h.angular_representation(1) + D3h.angular_representation(2) * D3h.irrep_characters[D3h.inversion_irrep]
print(f"Predicted coupled exciton rep: {D3h.reduce_representation(spherical_harmonic_prediction)[1]}")
print(f"Measured coupled exciton rep: {D3h.reduce_representation(D3h_X)[1]}")
print(f"Altered prediction assuming separability of geradeness: {D3h.reduce_representation(spherical_harmonic_ungerade)[1]}")


"""

The C3v -> D3h elevation study

both have rotoinversions, but C3v only has a 2-fold one (the three sigmas), so all the spinor reps have their characters be zero for it
Hence the multiplication of half-integer j angular reps by the inversion parity rep doesn't change anything

ONCE A GROUP CONTAINS A ROTOINVERSION WITH A FOLDEDNESS HIGHER THAN 2, THE SPINOR REPS ARE NO LONGER INVARIANT UNDER MULTIPLICATION BY THE PARITY REP.

basically not having that rotoinversion creates a degeneracy, as |1/2, +-1/2> and |3/2, +-1/2> form simultaneous bases.


"""
