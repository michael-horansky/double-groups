import matplotlib.pyplot as plt
import numpy as np


class ket():
    def __init__(self, degeneracy, indices):
        # the hilbert vector describing a D-fold degenerate energy level populated by a number of particles. len(indices) is the occupancy
        self.degeneracy = degeneracy
        self.i = indices
        self.j_z = [x - 1 - (degeneracy - 1) / 2 for x in self.i]
        self.occupancy = len(indices)
    
    def __repr__(self):
        if self.degeneracy % 2 == 1:
            return(f"|{','.join([str(int(a)) for a in self.j_z])}>")
        else:
            return(f"|{','.join([str(int(2*a))+'/2' for a in self.j_z])}>")
    
    def __str__(self):
        if self.degeneracy % 2 == 1:
            return(f"|{','.join([str(int(a)) for a in self.j_z])}>")
        else:
            return(f"|{','.join([str(int(2*a))+'/2' for a in self.j_z])}>")


def wigner_small_d_matrix(beta, j):
    dim = int(2*j+1)
    result = np.zeros((dim, dim), dtype=complex)
    for a in range(dim):
        m_a = -j + a
        for b in range(dim):
            m_b = -j + b
            partial_sum = 0.0
            for s in range(max(0, int(m_a - m_b)), min(int(j + m_a), int(j - m_b)) + 1):
                if m_b - m_a + s % 2 == 0:
                    pre_factor = 1
                else:
                    pre_factor = -1
                new_term = pre_factor * np.power(np.cos(beta / 2.0), 2*j+m_a-m_b-2*s) * np.power(np.sin(beta / 2.0), m_b-m_a+2*s)
                new_term /= (np.math.factorial(j + m_a - s) * np.math.factorial(s) * np.math.factorial(m_b - m_a + s) * np.math.factorial(j - m_b - s))
                partial_sum += new_term
            result[b][a] = partial_sum * np.sqrt(np.math.factorial(j+m_b) * np.math.factorial(j-m_b) * np.math.factorial(j+m_a) * np.math.factorial(j-m_a))
    return(result)

def direct_wigner_small_d_matrix(beta, j = 3/2):
    # this just directly implements d(beta) for j=3/2
    result = np.zeros((4, 4), dtype=complex)
    cos_b = np.cos(beta)
    cos_h_b = np.cos(beta / 2.0)
    sin_h_b = np.sin(beta / 2.0)
    sq_t = np.sqrt(3) / 2.0
    result[0][0] = 0.5 * (1 + cos_b) * cos_h_b
    result[0][1] = - sq_t * (1 + cos_b) * sin_h_b
    result[0][2] = sq_t * (1 - cos_b) * cos_h_b
    result[0][3] = -0.5 * (1 - cos_b) * sin_h_b
    
    result[1][0] = sq_t * (1 + cos_b) * sin_h_b
    result[1][1] = 0.5 * (3 * cos_b - 1) * cos_h_b
    result[1][2] = -0.5 * (3 * cos_b + 1) * sin_h_b
    result[1][3] = sq_t * (1 - cos_b) * cos_h_b
    
    result[2][0] = sq_t * (1 - cos_b) * cos_h_b
    result[2][1] = 0.5 * (3 * cos_b + 1) * sin_h_b
    result[2][2] = 0.5 * (3 * cos_b - 1) * cos_h_b
    result[2][3] = -sq_t * (1 + cos_b) * sin_h_b
    
    result[3][0] = 0.5 * (1 - cos_b) * sin_h_b
    result[3][1] = sq_t * (1 - cos_b) * cos_h_b
    result[3][2] = sq_t * (1 + cos_b) * sin_h_b
    result[3][3] = 0.5 * (1 + cos_b) * cos_h_b
    
    return(result)


def index_sublists(list_length, sublist_length, reverse_order = False):
    # returns a list of all lists of length sublist_length where elements are in strictly increasing order and for each 0<=x<list_length
    result = []
    if sublist_length == 1:
        for i in range(list_length):
            if reverse_order:
                result.append([list_length-i-1])
            else:
                result.append([i])
        return(result)
    if list_length == sublist_length:
        for i in range(list_length):
            if reverse_order:
                result.append(list_length-i-1)
            else:
                result.append(i)
        return([result])
    for i in range(list_length - sublist_length+1):
        partition_res = index_sublists(list_length - i - 1, sublist_length - 1)
        for part in partition_res:
            new_part = [i]
            for j in range(len(part)):
                if reverse_order:
                    new_part = [i + part[j] + 1] + new_part
                else:
                    new_part.append(i + part[j] + 1)
            result.append(new_part)
    return(result)


def find_representative_pure_states_by_occupancy(n, k):
    # find a single element from every permutation sublist
    all_index_sublists = index_sublists(n, k)
    representative_pure_states = []
    for i_sl in all_index_sublists:
        representative_pure_states.append(ket(n, [i + 1 for i in i_sl]))
    return(representative_pure_states)


def minor(matrix, indices):
    return(matrix[np.array(indices)[:,np.newaxis],np.array(indices)])

def get_occupancy_trace(k, theta, j):
    n = int(2*j+1)
    d_matrix = np.array(direct_wigner_small_d_matrix(theta, j))
    if k == 0:
        return(1)
    all_index_sublists = index_sublists(n, k)
    total_sum = 0.0
    for i_s in all_index_sublists:
        cur_d_minor = np.matrix(minor(d_matrix, i_s))
        total_sum += np.linalg.det(cur_d_minor)
    return(total_sum)
        

def get_all_occupancy_trace(theta, j):
    n = int(2*j+1)
    d_traces = [1]
    d_matrix = np.matrix(direct_wigner_small_d_matrix(theta, j))
    cur_matrix = np.matrix(d_matrix)
    for i in range(n):
        d_traces.append(np.trace(cur_matrix))
        cur_matrix = np.matmul(cur_matrix, d_matrix)
    results = [1]
    for m in range(1, n+1):
        cum_sum = 0.0
        for i in range(1, m+1):
            cum_sum += results[m-i] * d_traces[i]
        results.append(- cum_sum / m)
    return(results)

theta_space = np.linspace(0, np.pi, 100)

def rot_character(theta, j):
    return(np.sin((2*j+1) * theta / 2) / np.sin(theta / 2))

def actual_doublet(theta):
    return(2.0 * np.cos(theta) + 2.0 * np.cos(2.0 * theta) + 2)

a = [[], [], [], [], []] #i-th element is the trace of occupancy i
for t in theta_space:
    res = get_all_occupancy_trace(t, 3/2)
    for i in range(len(res)):
        a[i].append(res[i])
b = [[], [], [], [], []] #i-th element is the trace of occupancy i
for t in theta_space:
    for i in range(4+1):
        res = get_occupancy_trace(i, t, 3/2)
        b[i].append(res)
        

lol = find_representative_pure_states_by_occupancy(4, 2)
print(lol)

plt.xlabel("theta")
plt.ylabel("character")
for i in range(len(a)):
    plt.plot(theta_space, b[i], label=f"occupancy {i}")
plt.plot(theta_space, actual_doublet(theta_space), linestyle="dashed", label="actual doublet")
plt.plot(theta_space, rot_character(theta_space, 0)+rot_character(theta_space, 2), linestyle="dashed", label="decomposed doublet")
plt.legend()
plt.show()

#def find_antisymmetric_basis(n, k):
#    # finds the antisymmetric basis of an n-fold degenerate energy level with occupancy k, where the i-th basis vector is labelled as |i>, i = 1 ..., n-1

