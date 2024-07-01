import matplotlib.pyplot as plt
import numpy as np
import scipy as sp


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
            cum_sum += (1-(i % 2) * 2) * results[m-i] * d_traces[i]
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
        

def assumed_wigner_power_trace(theta, k):
    j = 3/2
    return(np.sin((2 * j + 1) * k * theta / 2) / np.sin(k * theta / 2))

"""k_space = np.arange(1, 10)
cur_theta = 1.0
actual_result_space = []
cur_d_matrix = np.matrix(direct_wigner_small_d_matrix(cur_theta, 3/2))
cur_d_matrix_e = np.matrix(cur_d_matrix)
for k_val in k_space:
    actual_result_space.append(np.trace(cur_d_matrix_e))
    cur_d_matrix_e = np.matmul(cur_d_matrix_e, cur_d_matrix)

plt.plot(k_space, actual_result_space, label = "actual d^k")
plt.plot(k_space, assumed_wigner_power_trace(cur_theta, k_space), linestyle="dashed", label = "assumed d^k")
plt.show()"""
    

"""lol = find_representative_pure_states_by_occupancy(4, 2)
print(lol)

plt.xlabel("theta")
plt.ylabel("character")
for i in range(len(a)):
    plt.plot(theta_space, a[i], label=f"occupancy {i}")
plt.plot(theta_space, actual_doublet(theta_space), linestyle="dashed", label="actual doublet")
plt.plot(theta_space, rot_character(theta_space, 0)+rot_character(theta_space, 2), linestyle="dashed", label="decomposed doublet")
plt.legend()
plt.show()"""

#def find_antisymmetric_basis(n, k):
#    # finds the antisymmetric basis of an n-fold degenerate energy level with occupancy k, where the i-th basis vector is labelled as |i>, i = 1 ..., n-1



def bitfield(n, l):
    #return [1 if digit=='1' else 0 for digit in bin(n)[2:]]
    return [1 if digit=='1' else 0 for digit in  f'{n:0{l}b}']
def bitfield_to_composition(bf):
    res = []
    cur_sum = 1
    for i in range(len(bf)):
        if bf[i] == 1:
            cur_sum += 1
        else:
            res.append(cur_sum)
            cur_sum = 1
    #if cur_sum > 1:
    res.append(cur_sum)
    return(res)

def compositions_of_k(k):
    if k == 0:
        return([])
    if k == 1:
        return([[1]])
    res = []
    for i in range(np.power(2, k-1)):
        cur_bf = bitfield_to_composition(bitfield(i, k-1))
        res.append(cur_bf)
    return(res)


def get_Omega(b, j):
    k = sum(b)
    # the lower non-zero value of Omega is at -j * sum(b) = -j*k, with the upper bound at j*k. Hence there are 4jk+1 possible values of omega for which Omega is non-zero. These will be returned as a list with 4jk+1 elements, indices ranging from 0 to 4jk, with the property that Omega[2jk] = Omega(0), and Omega[4jk-i]=Omega[i].
    i_null = int(2 * j * k)
    Omega = [0] * (i_null * 2+1)
    def omega_from_index(i):
        return((i - i_null) / 2)
    def index_from_omega(omega):
        return(int(2 * omega) + i_null)
    
    def add_next_element_to_Omega(cur_Omega, new_index):
        next_Omega = [0] * len(cur_Omega)
        for i in range(len(cur_Omega)):
            if cur_Omega[i] > 0:
                cur_omega = omega_from_index(i)
                for c in np.arange(-j, j+1, 1):
                    new_omega = cur_omega + c * b[new_index]
                    next_Omega[index_from_omega(new_omega)] += cur_Omega[i]
        return(next_Omega)
    
    # we start with 1 and then iterate over every sum generated by b_i
    Omega[i_null] = 1
    for i in range(len(b)):
        #print(Omega)
        Omega = add_next_element_to_Omega(Omega, i)
    return(Omega)

def decomposition_of_multifermionic_ensemble(j, k):
    # returns a list c_jk[i], where c_jk[2j'] is the coef of c^j_k(j').
    # we choose the upper bound on j' as (degeneracy-1)/2, where degeneracy is the number of antisymmetric basis vectors
    degeneracy = sp.special.binom(int(2 * j +1), k)
    max_j = (degeneracy - 1)/2
    i_max = int(max_j * 2)
    
    #c_jk = [0] * (i_max + 1)
    c_jk = [0] * (i_max + 1)*k
    
    def sign_exp(x):
        if x % 2 == 0:
            return(1)
        else:
            return(-1)
    
    k_comps = compositions_of_k(k)
    
    for k_c in k_comps:
        #print(k_c)
        n = len(k_c)
        # first, we calculate the product factor
        product_factor = k
        for i in range(1, n):
            product_factor *= (k - sum(k_c[:i]))
        # now, we calculate Omega
        Omega = get_Omega(k_c, j)
        #print(Omega)
        #print("Product factor = 1 /", product_factor)
        # then, we do the thing!
        for i in range(int(2 * j * k), len(Omega)-2):
            cur_j = (i - 2 * j * k) / 2
            #if cur_j > max_j:
            #    continue
            c_jk[int(2 * cur_j)] += sign_exp(k+n) / product_factor * (Omega[i] - Omega[i+2])
        for i in range(len(Omega)-2, len(Omega)):
            cur_j = (i - 2 * j * k) / 2
            #if cur_j > max_j:
            #    continue
            #print(cur_j)
            c_jk[int(2 * cur_j)] += sign_exp(k+n) / product_factor * (Omega[i])
    
    # now we round
    for i in range(len(c_jk)):
        c_jk[i] = int(np.round(c_jk[i]))
    
    return(c_jk)



# ------------------- the recursive method ----------------------

def Decomposition_from_product(x, y):
    new_keys = np.arange(np.abs(x-y), x+y+1, 1)
    new_dict = {}
    for k in new_keys:
        new_dict[k] = 1
    return(Decomposition(new_dict))

class Decomposition():
    
    def __init__(self, nonzero_coefs):
        # nonzero_coefs is a dictionary = {float j : int c_j}
        self.nonzero_coefs = nonzero_coefs
    
    def __str__(self):
        return(str(self.nonzero_coefs))
    
    def coef(self, j):
        if j in self.nonzero_coefs.keys():
            return(int(self.nonzero_coefs[j]))
        else:
            return(0)
    
    def add_to_coef(self, j, x):
        # we increase self.coef(j) by x
        if j in self.nonzero_coefs.keys():
            self.nonzero_coefs[j] += x
        else:
            self.nonzero_coefs[j] = x
        
    def __add__(self, other):
        new_nonzero_coefs = {}
        for my_key in self.nonzero_coefs.keys():
            new_nonzero_coefs[my_key] = self.nonzero_coefs[my_key]
        for other_key in other.nonzero_coefs.keys():
            if other_key in new_nonzero_coefs.keys():
                new_nonzero_coefs[other_key] += other.nonzero_coefs[other_key]
            else:
                new_nonzero_coefs[other_key] = other.nonzero_coefs[other_key]
        return(Decomposition(new_nonzero_coefs))
    
    def __mul__(self, other):
        final_result = Decomposition({})
        for my_key in self.nonzero_coefs.keys():
            for other_key in other.nonzero_coefs.keys():
                final_result += Decomposition_from_product(my_key, other_key).scalar_mul(self.coef(my_key) * other.coef(other_key))
        return(final_result)
    
    def scalar_mul(self, c):
        result = Decomposition(self.nonzero_coefs.copy())
        for k in result.nonzero_coefs.keys():
            result.nonzero_coefs[k] *= c
        return(result)
    
    def mul_with_Dirichlet_contraction(self, a):
        # multiply self with chi^j - chi^{j-1}
        
        # the A to B method
        result = Decomposition({})
        for b in np.arange(0, a-1/2, 1/2):
            result.add_to_coef(a+b, self.coef(b))
            result.add_to_coef(a-b-1, -self.coef(b))
        result.add_to_coef(2*a-1/2, self.coef(a-1/2))
        result.add_to_coef(0, self.coef(a))
        result.add_to_coef(2*a, self.coef(a))
        for b in np.arange(a+1/2, max(self.nonzero_coefs.keys())+1/2, 1/2):
            result.add_to_coef(a+b, self.coef(b))
            result.add_to_coef(b-a, self.coef(b))
        return(result)


def decomposition_of_multifermionic_ensemble_recursively(j, k):
    # returns a decomposition
    
    # edge cases
    if j == 0:
        return(Decomposition({0 : 1}))
    if k == 0 or k == int(2 * j + 1):
        return(Decomposition({0 : 1}))
    if k == 1 or k == int(2 * j):
        return(Decomposition({j : 1}))
    
    # tail recursion
    left_part  = decomposition_of_multifermionic_ensemble_recursively(j - 1/2, k  ).mul_with_Dirichlet_contraction(k/2      )
    right_part = decomposition_of_multifermionic_ensemble_recursively(j - 1/2, k-1).mul_with_Dirichlet_contraction((2*j+1-k)/2)
    return((left_part + right_part).scalar_mul(1/2))






        
def print_table(table_name, column_names, row_names, list_of_rows, subtable_borders = [], header_separation = 2, empty_cell_str=" ", omit_strings=[], print_to_stdout = True):
    # column_names[N], row_names[M], list_of_rows[M][N]
    # for each element in row_names which is just "-", a horizontal sub-border is printed instead
    def st(a, w):
        if len(str(a)) >= w:
            return(str(a))
        else:
            diff = w - len(str(a))
            return(" " * int(np.floor(diff / 2.0)) + str(a) + " " * int(np.ceil(diff / 2.0)))
    
    result_string = ""
    
    max_len_row_names = len(str(table_name))
    for row_name in row_names:
        if len(str(row_name)) > max_len_row_names:
            max_len_row_names = len(row_name)
    max_len_by_column = []
    for column_name in column_names:
        max_len_by_column.append(len(str(column_name)))
    
    for j in range(len(column_names)):
        skipped_lines = 0
        for i in range(len(row_names)):
            if row_names[i] == "-":
                skipped_lines += 1
                continue
            if len(list_of_rows[i-skipped_lines]) > j:
                if max_len_by_column[j] < len(str(list_of_rows[i-skipped_lines][j])):
                    max_len_by_column[j] = len(str(list_of_rows[i-skipped_lines][j]))
    
    header_str = st(table_name, max_len_row_names + header_separation) + "|"
    for i in range(len(column_names)):
        header_str += st(column_names[i], max_len_by_column[i] + header_separation)
        if i in subtable_borders:
            header_str += "|"
    
    skipped_lines = 0
    
    result_string = header_str + "\n" + "-" * len(header_str) + "\n"
    for i in range(len(row_names)):
        if row_names[i] == "-":
            skipped_lines += 1
            result_string += "-" * len(header_str) + "\n"
            continue
        cur_str = st(row_names[i], max_len_row_names + header_separation) + "|"
        for j in range(len(column_names)):
            if len(list_of_rows[i-skipped_lines]) > j:
                if str(list_of_rows[i-skipped_lines][j]) in omit_strings:
                    cur_str += st("", max_len_by_column[j] + header_separation)
                else:
                    cur_str += st(list_of_rows[i-skipped_lines][j], max_len_by_column[j] + header_separation)
            else:
                cur_str += st(empty_cell_str, max_len_by_column[j] + header_separation)
            if j in subtable_borders:
                cur_str += "|"
        result_string += cur_str + "\n"
    result_string += "-" * len(header_str) + "\n"
    if print_to_stdout:
        print(result_string)
    return(result_string)




#cur_b = [3, 1, 1] #NOTE these are just convolutions with itself!! discrete convolution? dirichlet kernel properties?
#cur_j = 3/2

#print(compositions_of_k(1))

#print(get_Omega(cur_b, cur_j))

#print("viola", decomposition_of_multifermionic_ensemble(3/2, 3))






#cur_j = 3/2
#max_row_length = int(sp.special.binom(int(2 * 2 +1), int(cur_j+1/2)))

def decomposition_table_column_names(row_length):
    res = []
    for i in range(row_length):
        if i % 2 == 0:
            res.append(str(int(i/2)))
        else:
            res.append(str(i) + "/2")
    return(res)

"""decomposition_table = []
decomposition_table_row_names = []

for k in range(1, int(cur_j * 2 + 1)):
    cur_decomposition = decomposition_of_multifermionic_ensemble(cur_j, k)
    decomposition_table.append(cur_decomposition[:max_row_length])
    if int(2*cur_j) % 2 == 0:
        decomposition_table_row_names.append(str(int(cur_j)) + ", " + str(k))
    else:
        decomposition_table_row_names.append(str(int(2*cur_j)) + "/2, " + str(k))
cur_j = 2

decomposition_table_row_names.append("-")

for k in range(1, int(cur_j * 2 + 1)):
    cur_decomposition = decomposition_of_multifermionic_ensemble(cur_j, k)
    decomposition_table.append(cur_decomposition[:max_row_length])
    if int(2*cur_j) % 2 == 0:
        decomposition_table_row_names.append(str(int(cur_j)) + ", " + str(k))
    else:
        decomposition_table_row_names.append(str(int(2*cur_j)) + "/2, " + str(k))"""


#print_table("c^j_k", decomposition_table_column_names(max_row_length), decomposition_table_row_names, decomposition_table, omit_strings=["0"])

def print_decomposition_table(list_of_j):
    
    decomposition_table_row_names = []
    decomposition_table = []
    
    actual_max_length = 0
    
    for cur_j in list_of_j:
        print("Working on j =", cur_j)
        for k in range(1, int(cur_j * 2 + 1)):
            cur_decomposition = decomposition_of_multifermionic_ensemble(cur_j, k)
            
            # we find the position of the largest nonzero element
            for i in range(len(cur_decomposition)):
                if cur_decomposition[len(cur_decomposition)-1-i] != 0:
                    break
            cur_max_length = len(cur_decomposition)-i
            if actual_max_length < cur_max_length:
                actual_max_length = cur_max_length
            
            decomposition_table.append(cur_decomposition[:cur_max_length])
            if int(2*cur_j) % 2 == 0:
                decomposition_table_row_names.append(str(int(cur_j)) + ", " + str(k))
            else:
                decomposition_table_row_names.append(str(int(2*cur_j)) + "/2, " + str(k))
        #print("    New max length =", actual_max_length)
        if cur_j != list_of_j[-1]:
            decomposition_table_row_names.append("-")
    
    # we trim everything
    for i in range(len(decomposition_table)):
        decomposition_table[i] = decomposition_table[i][:actual_max_length]
    
    print_table("c^j_k", decomposition_table_column_names(actual_max_length), decomposition_table_row_names, decomposition_table, omit_strings=["0"])

#print_decomposition_table(np.arange(1/2, 8/2, 1/2))

def print_decomposition_table_from_recursion(list_of_j):
    
    decomposition_table_row_names = []
    decomposition_table = []
    
    actual_max_length = 0
    
    for cur_j in list_of_j:
        print("Working on j =", cur_j)
        for k in range(1, int(cur_j * 2 + 1)):
            cur_decomposition = decomposition_of_multifermionic_ensemble_recursively(cur_j, k)
            
            max_j = max(cur_decomposition.nonzero_coefs.keys())
            
            cur_tab_row = []
            for j in np.arange(0, max_j+1/2, 1/2):
                cur_tab_row.append(cur_decomposition.coef(j))
            cur_max_length = len(cur_tab_row)
            if actual_max_length < cur_max_length:
                actual_max_length = cur_max_length
            
            decomposition_table.append(cur_tab_row)
            if int(2*cur_j) % 2 == 0:
                decomposition_table_row_names.append(str(int(cur_j)) + ", " + str(k))
            else:
                decomposition_table_row_names.append(str(int(2*cur_j)) + "/2, " + str(k))
        #print("    New max length =", actual_max_length)
        if cur_j != list_of_j[-1]:
            decomposition_table_row_names.append("-")
    
    # we trim everything
    for i in range(len(decomposition_table)):
        decomposition_table[i] = decomposition_table[i][:actual_max_length]
    
    print_table("c^j_k", decomposition_table_column_names(actual_max_length), decomposition_table_row_names, decomposition_table, omit_strings=["0"])


print_decomposition_table_from_recursion(np.arange(1/2, 12/2, 1/2)) # a massive time improvement, but we should be saving the tail recursion results, which axes the master complexity from square to linear!
