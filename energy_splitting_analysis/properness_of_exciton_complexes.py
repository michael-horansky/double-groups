


from aretdog.class_QD_group import *


# Here we will test the differences between the parity-propagated and the naive j-irreps and their allowed transitions

C3v = QDGroup("C3v QD")
C3v.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "sigma" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True)})
C3v.print_character_table()
C3v.classify_excitons(max_e = 2, max_h = 2)
C3v.print_exciton_complexes()

C3v.find_transition_chain()
print(C3v.transition_chain)

D3h = QDGroup("D3h QD")
D3h.generate_double_group({"E" : E, "Cz_3" : ImproperRotation([0.0, 0.0, 1.0], [1, 3], False), "m" : ImproperRotation([1.0, 1.0, 0.0], [1, 2], True), "m'" : ImproperRotation([0.0, 0.0, 1.0], [1, 2], True)})
D3h.print_character_table()
D3h.classify_excitons(max_e = 2, max_h = 2)
D3h.print_exciton_complexes()

D3h.find_transition_chain()
print(D3h.transition_chain)


print("----------------- Predictions of spectral features -------------------")

# e: j=1/2; h: j=3/2. Hence X_10+X_01 = Gamma_1 + Gamma_2
# 2X_11: e+e = identity; h1 + h2 = h1 x h2
print("------------ C3v group")
#print("natural j=1: ", C3v.angular_representation(1))
#print("natural j=2: ", C3v.reduce_representation(C3v.angular_representation(2))[1])

print("----Naive representation")

#C3v_X_1_naive = C3v.angular_representation(1) + C3v.angular_representation(2, "u")
C3v_X_1_naive = C3v.angular_representation(1) + C3v.angular_representation(2)
C3v_2X_11_naive = C3v.excitons["1h1"] * C3v.excitons["1h2"]
print("X_10 + X_01 = ", C3v.reduce_representation(C3v_X_1_naive)[1])
print("2X_11 = ", C3v.reduce_representation(C3v_2X_11_naive)[1])

C3v_2X_transitions_naive = C3v.allowed_transitions_between_reps(C3v_2X_11_naive, C3v_X_1_naive)
C3v_ground_transitions_naive = C3v.allowed_transitions_between_reps(C3v_X_1_naive, C3v.excitons["vacuum"])

print("Allowed transitions (2X_11 -> X_10) + (2X_01 -> X_01):")

total_transitions = 0
for polaris, trans in C3v_2X_transitions_naive.items():
    print(f"  {polaris}-polar.: {'; '.join(trans[0])}")
    total_transitions += len(trans[0])
print("  Total allowed transitions:", total_transitions)

print("Allowed transitions (X_10 -> vacuum) + (X_01 -> vacuum):")

total_transitions = 0
for polaris, trans in C3v_ground_transitions_naive.items():
    print(f"  {polaris}-polar.: {'; '.join(trans[0])}")
    total_transitions += len(trans[0])
print("  Total allowed transitions:", total_transitions)

print("")
print("----Parity-conscious representation")

C3v_X_1 = C3v.excitons["1X[1][1,0]"] + C3v.excitons["1X[1][0,1]"]
C3v_2X_11 = C3v.excitons["2X[2][1,1]"]
print("X_10 + X_01 = ", C3v.reduce_representation(C3v_X_1)[1])
print("2X_11 = ", C3v.reduce_representation(C3v_2X_11)[1])

C3v_2X_transitions = C3v.allowed_transitions_between_reps(C3v_2X_11, C3v_X_1)
C3v_ground_transitions = C3v.allowed_transitions_between_reps(C3v_X_1, C3v.excitons["vacuum"])


print("Allowed transitions (2X_11 -> X_10) + (2X_01 -> X_01):")

total_transitions = 0
for polaris, trans in C3v_2X_transitions.items():
    print(f"  {polaris}-polar.: {'; '.join(trans[0])}")
    total_transitions += len(trans[0])
print("  Total allowed transitions:", total_transitions)

print("Allowed transitions (X_10 -> vacuum) + (X_01 -> vacuum):")

total_transitions = 0
for polaris, trans in C3v_ground_transitions.items():
    print(f"  {polaris}-polar.: {'; '.join(trans[0])}")
    total_transitions += len(trans[0])
print("  Total allowed transitions:", total_transitions)
print("")
"""
print("------------ D3h group")

print("----Naive representation")

D3h_X_1_naive = D3h.angular_representation(1) + D3h.angular_representation(2)
D3h_2X_11_naive = D3h.excitons["1h1"] * D3h.excitons["1h2"]
print("X_10 + X_01 = ", D3h.reduce_representation(D3h_X_1_naive)[1])
print("X_11 = ", D3h.reduce_representation(D3h_2X_11_naive)[1])

D3h_2X_transitions_naive = D3h.allowed_transitions_between_reps(D3h_2X_11_naive, D3h_X_1_naive)
D3h_ground_transitions_naive = D3h.allowed_transitions_between_reps(D3h_X_1_naive, D3h.excitons["vacuum"])

print("Allowed transitions (2X_11 -> X_10) + (2X_01 -> X_01):")

total_transitions = 0
for polaris, trans in D3h_2X_transitions_naive.items():
    print(f"  {polaris}-polar.: {'; '.join(trans[0])}")
    total_transitions += len(trans[0])
print("  Total allowed transitions:", total_transitions)
print("Allowed transitions (X_10 -> vacuum) + (X_01 -> vacuum):")

total_transitions = 0
for polaris, trans in D3h_ground_transitions_naive.items():
    print(f"  {polaris}-polar.: {'; '.join(trans[0])}")
    total_transitions += len(trans[0])
print("  Total allowed transitions:", total_transitions)

print("")
print("----Parity-conscious representation")

D3h_X_1 = D3h.excitons["1X[1][1,0]"] + D3h.excitons["1X[1][0,1]"]
D3h_2X_11 = D3h.excitons["2X[2][1,1]"]
print("X_10 + X_01 = ", D3h.reduce_representation(D3h_X_1)[1])
print("X_11 = ", D3h.reduce_representation(D3h_2X_11)[1])

D3h_2X_transitions = D3h.allowed_transitions_between_reps(D3h_2X_11, D3h_X_1)
D3h_ground_transitions = D3h.allowed_transitions_between_reps(D3h_X_1, D3h.excitons["vacuum"])


print("Allowed transitions (2X_11 -> X_10) + (2X_01 -> X_01):")

total_transitions = 0
for polaris, trans in D3h_2X_transitions.items():
    print(f"  {polaris}-polar.: {'; '.join(trans[0])}")
    total_transitions += len(trans[0])
print("  Total allowed transitions:", total_transitions)

print("Allowed transitions (X_10 -> vacuum) + (X_01 -> vacuum):")

total_transitions = 0
for polaris, trans in D3h_ground_transitions.items():
    print(f"  {polaris}-polar.: {'; '.join(trans[0])}")
    total_transitions += len(trans[0])
print("  Total allowed transitions:", total_transitions)
"""
