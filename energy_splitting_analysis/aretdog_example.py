


from aretdog.class_group import *

C_i = Group("C_i")
C_i.generate_group({"E" : E, "i" : ImproperRotation([1.0, 0.0, 0.0], [0, 1], True)})
C_i.print_character_table()


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


"""
Th_group = Group("Th")
Th_group.generate_double_group({"E" : E, "Cz_2" : Cz_2, "Cy_2" : Cy_2, "C'_3" : ImproperRotation([1.0, 1.0, 1.0], [1, 3], False), "i" : ImproperRotation([1.0, 0.0, 0.0], [0.0, 1.0], True)})
Th_group.print_character_table()"""



print("------------------ C3v analysis ------------------------")

C3v_group = Group("C3v")
C3v_group.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C3v_group.print_character_table()




D3h_group = Group("D3h")
D3h_group.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "m" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True), "m'" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], True)})
D3h_group.print_character_table()


D3h_group.add_subgroup(C3v_group)

two_X_11 = C3v_group.irrep_characters["E_1(j=1/2)"] * (C3v_group.irrep_characters["1E(j=3/2)"] + C3v_group.irrep_characters["2E(j=3/2)"])
X_10 = C3v_group.irrep_characters["E_1(j=1/2)"] * (C3v_group.irrep_characters["1E(j=3/2)"] + C3v_group.irrep_characters["2E(j=3/2)"])
X_01 = C3v_group.irrep_characters["E_1(j=1/2)"] * C3v_group.irrep_characters["E_1(j=1/2)"]
vacuum = C3v_group.irrep_characters["A_1"]
C3v_dipole_rep = C3v_group.irrep_characters["E_1"] # Karlsson p 13 - in an optic transition, we do dipole approximation and Wigner-Eckart (at least for isotropic polarisation)

print("2X_11 has energy leveles as", C3v_group.reduce_representation(two_X_11)[1])
print("X_10 has energy leveles as", C3v_group.reduce_representation(X_10)[1])
print("X_01 has energy leveles as", C3v_group.reduce_representation(X_01)[1])

X_11_to_X_10 = C3v_group.allowed_transitions_between_reps(two_X_11, X_10)

for polarisation in X_11_to_X_10:
    print(polarisation + "-polarised transitions in 2X_11 -> X_10:")
    print("  Allowed: ", X_11_to_X_10[polarisation][0])
    print("  Dark:    ", X_11_to_X_10[polarisation][1])

X_01_to_vacuum = C3v_group.allowed_transitions_between_reps(X_01, vacuum)
for polarisation in X_01_to_vacuum:
    print(polarisation + "-polarised transitions in X_01 -> vacuum:")
    print("  Allowed: ", X_01_to_vacuum[polarisation][0])
    print("  Dark:    ", X_01_to_vacuum[polarisation][1])



print("------------------ exciton complex labelling analysis ------------------------")


print("C3v:")




print("D3h:")

# if hole = (spinless wavefunction) * (j=1/2 wavefunction) = E_j=5/2, what is that spinless wavefunction?

print("s-orbital irrep :", D3h_group.reduce_representation(D3h_group.angular_representation(0))[1])
print("p-orbital irrep :", D3h_group.reduce_representation(D3h_group.angular_representation(1))[1])

print("p-orbital x j=3/2:", D3h_group.reduce_representation(D3h_group.angular_representation(1) * D3h_group.angular_representation(3/2))[1])

print("p-orbital x j=3/2 MINUS j=3/2 (split-off):", D3h_group.reduce_representation(D3h_group.angular_representation(1) * D3h_group.angular_representation(3/2) - D3h_group.angular_representation(3/2))[1])

# this is bullshit since the reason hole has j=3/2 is because we coupled l=1 with s=1/2. This is double coupling = wrong
print("p-orbital x j=1/2:", D3h_group.reduce_representation(D3h_group.angular_representation(1) * D3h_group.angular_representation(1/2))[1])



print("A_3 * E(j=1/2) = ", D3h_group.reduce_representation(D3h_group.irrep_characters["A_3"]*D3h_group.irrep_characters["E_2(j=1/2)"])[1])
print("A_3 * E(j=3/2) = ", D3h_group.reduce_representation(D3h_group.irrep_characters["A_3"]*D3h_group.irrep_characters["E_1(j=3/2)"])[1])

print("Hence ASSUMING the envelope symmetry of holes is always the inversion parity rep, we could construct hole states like", D3h_group.reduce_representation(D3h_group.irrep_characters["A_3"]*D3h_group.angular_representation(3/2))[1])

# THIS IS IT!! the spin parts of the hole wavefunctions are the j=3/2, j_z = (+-3/2) for HH and (+-1/2) for LH
# however there's also the spatial envelope function A_3 (A_1'' in Altmann) which we gotta multiply that with!!
# But then... why does the electron spatial envelope func transform according to Identity, not A_3? Hmmm thonk
# remember electrons are s orbital and holes are p orbital!!!!!

# p-orbital transforms according to A_4 + E_1, NOT A_3.... how?


# Notes on exciton complex representations:
"""

Karlsson p. 5: "There is no spectral evidence of any additional hole level confined to the QD, besides the discussed hh- and
lh-like levels. Moreover, there are no indications of any excited electron levels in the spectra. Consequently, all
the details of the pyramidal QD spectra will in the following be analyzed on the basis of one single electron level,
one hh-like level and one lh-like level. Each level is assumed to accommodate two degenerate single particle
states due to Kramers degeneracy."

This basically means that electrons are assumed to have l=0 -> j=l+s=1/2, which why they always transform as E_1/2


Dalessi, p. 10: "For the valence band of semiconductors
like GaAs, at the so-called GAMMA point, one restricts to j = 3 / 2
(Ref. 22) giving rise to the Luttinger Hamiltonian [cf. Eqs.
(3) and (4)]. However, the formalism holds for an arbitrary
large number of bands."

remember a hole is just a missing electron state at the top of the valence band!! this is it!!!


Dalessi further states that HH are associated with Jz=+-3/2 and LH are associated with Jz=+-1/2. How do we find what these transform like? Easy!
Take the full J=3/2 angular rep. This rep must have the four spherical harmonic Y_lm, l=3/2 as its bases. Proceed analogously with (x,y,z) to see which
irreps transform like what spherical harmonics. NOTE you might need to consider "physical" irreps only (combine complex onjugate irrep pairs by summing)


"""





