
import numpy as np


def numpy_gcd(a, b):
    while b != 0:
        t = b
        b = a % b
        a = t
    return a

def reduce_frac(frac):
    if frac[0] == 0:
        return(frac)
    gcd = numpy_gcd(frac[0], frac[1])
    return(constrained_angle([int(frac[0] / gcd), int(frac[1] / gcd)]))

def m_to_angle(multiplicity):
    return(2.0 * np.pi * multiplicity[0] / multiplicity[1])

def float_to_fraction(x, max_q = 100):
    best_q = -1
    best_p = -1
    best_offset = 1e9
    for q in range(1, max_q + 1):
        cur_numer = np.round(x * q)
        val = cur_numer/q
        offset = abs(1 - val/x)
        if offset < best_offset:
            best_q = q
            best_p = int(cur_numer)
            best_offset = offset
    return([best_p, best_q])


def angle_to_m(angle, max_q = 12):
    # nontrivial - a rational fraction approximator
    
    while angle < 0.0:
        angle += 2.0 * np.pi
    while angle >= 2.0 * np.pi:
        angle -= 2.0 * np.pi
    
    if np.round(angle, 10) == 0.0:
        return([0, 1])
    
    x = angle / (2.0 * np.pi)
    
    
    return(constrained_angle(reduce_frac(float_to_fraction(x, max_q))))
    

def constrained_angle(multiplicity):
    # constrains p into the [0, q-1] interval
    res = multiplicity[0]
    while(res < 0):
        res += multiplicity[1]
    while(res >= multiplicity[1]):
        res -= multiplicity[1]
    return([res, multiplicity[1]])

def inverse_angle(multiplicity):
    # given an angle around an axis, this returns the positive angle within (0, 2pi) which inverts that rotation
    res = multiplicity.copy()
    if res[0] != 0:
        res[0] = res[1] - res[0]
    return(constrained_angle(res))

def add_multiplicities(m1, m2):
    new_p = m1[0] * m2[1] + m2[0] * m1[1]
    new_q = m1[1] * m2[1]
    return(reduce_frac([new_p, new_q]))

def seitz_notation(axis, foldedness, inversion):
    
    sign_list = []
    rational_axis = np.array(axis)
    for i in range(len(rational_axis)):
        sign_list.append(rational_axis[i] < 0.0)
        if rational_axis[i] == 0.0:
            continue
        rational_axis /= rational_axis[i]
    
    list_of_fractions = [[], []]
    for i in range(len(rational_axis)):
        if rational_axis[i] == 0.0:
            continue
        cur_frac = float_to_fraction(rational_axis[i], 12)
        list_of_fractions[0].append(cur_frac[0])
        list_of_fractions[1].append(cur_frac[1])
    q_gcd = np.gcd.reduce(list_of_fractions[1])
    p_gcd = np.gcd.reduce(list_of_fractions[0])
    product_factor = np.prod(list_of_fractions[1])
    factor = product_factor / (q_gcd * p_gcd)
    str_axis = ""
    j = 0
    for i in range(len(rational_axis)):
        if sign_list[i]:
            str_axis += "-"
        if rational_axis[i] == 0.0:
            str_axis += "0"
        else:
            str_axis += str(int(factor * list_of_fractions[0][j] / list_of_fractions[1][j]))
            j += 1
    
    str_operation = ""
    
    if inversion:
        # either a reflection, an inversion, or a rotoinversion
        if foldedness == 2:
            # reflextion
            str_operation = "m"
        else:
            str_operation = "-" + str(int(foldedness))
    else:
        str_operation = str(int(foldedness))
    return(f"{str_operation}_{str_axis}")


def Euler_Rodriguez_encoder(axis, angle):
    a = np.cos(angle / 2.0)
    b = axis[0] * np.sin(angle / 2.0)
    c = axis[1] * np.sin(angle / 2.0)
    d = axis[2] * np.sin(angle / 2.0)
    return(a, b, c, d)

def Euler_Rodriguez_decoder(a, b, c, d):
    axis = np.zeros(3)
    angle = 2.0 * np.arccos(a)
    while(angle < 0.0):
        angle += 2.0 * np.pi
    while(angle >= 2.0 * np.pi):
        angle -= 2.0 * np.pi
    
    # special case: for a null angle, the axis is arbitrary
    if angle == 0.0:
        axis[0] = 1.0
        return(axis, angle)
    
    sin_angle_half = np.sqrt(1.0 - a * a)
    axis[0] = b / sin_angle_half
    axis[1] = c / sin_angle_half
    axis[2] = d / sin_angle_half
    return(axis, angle)




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


def integer_partitions(total_sum, number_of_constituents):
    # returns a list where each element is a list of length number_of_constituents and sum total_sum
    # fermionic occupancies for a given macrostate
    if total_sum == 0:
        return([[0] * number_of_constituents])
    if number_of_constituents == 1:
        return([[total_sum]])
    else:
        result = []
        for leftmost_element in range(0, total_sum + 1):
            cur_partitions = integer_partitions(total_sum - leftmost_element, number_of_constituents - 1)
            for partition in cur_partitions:
                result.append([leftmost_element] + partition)
        return(result)


