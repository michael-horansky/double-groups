import matplotlib.pyplot as plt
import numpy as np



# ------------------- printing functions ---------------------
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

def decomposition_table_column_names(row_length):
    res = []
    for i in range(row_length):
        if i % 2 == 0:
            res.append(str(int(i/2)))
        else:
            res.append(str(i) + "/2")
    return(res)



# -------------------------- calculation ---------------------------

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
    
    def scalar_mul(self, c):
        result = Decomposition(self.nonzero_coefs.copy())
        for k in result.nonzero_coefs.keys():
            result.nonzero_coefs[k] *= c
        return(result)
        
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

def print_decomposition_table_from_recursion(j_upper_cap, test_recursive_relation = False):
    
    decomposition_table_row_names = []
    decomposition_table = []
    
    actual_max_length = 0
    
    master_list = [] # master_list[j_index][k] = [c_0, c_1/2, c_1... c_max_j], where j_index = j * 2, so the first item corresponds to j = 0, and k = 0... 2j+1
    for j_index in range(int(j_upper_cap * 2)+1):
        
        cur_j = j_index / 2
        master_list.append([0] * (j_index + 2))
        master_list[j_index][0] = Decomposition({0 : 1})
        master_list[j_index][j_index + 1] = master_list[j_index][0]
        master_list[j_index][1] = Decomposition({cur_j : 1})
        master_list[j_index][j_index] = master_list[j_index][1]
        for k in range(2, int(np.ceil(cur_j + 1))):
            left_part  = master_list[j_index-1][k  ].mul_with_Dirichlet_contraction(k/2      )
            right_part = master_list[j_index-1][k-1].mul_with_Dirichlet_contraction((j_index+1-k)/2)
            master_list[j_index][k] = (left_part + right_part).scalar_mul(1/2)
            master_list[j_index][j_index + 1 - k] = master_list[j_index][k]
            
            if test_recursive_relation:
                # NOTE here we test the formulation from the article, which is inverse to the one we use here (so we test the validity of the inversion)
                max_j = max(master_list[j_index][k].nonzero_coefs.keys())
                
                for j_prime in np.arange(0, max_j+1/2, 1/2):
                    cur_coef = master_list[j_index][k].coef(j_prime)
                    # alpha
                    if j_prime <= k/2 - 1:
                        alpha_j_prime = master_list[j_index-1][k].coef(k/2 + j_prime) - master_list[j_index-1][k].coef(k/2 - j_prime - 1)
                    elif j_prime == (k-1)/2:
                        alpha_j_prime = master_list[j_index-1][k].coef(k - 1/2)
                    else:
                        alpha_j_prime = master_list[j_index-1][k].coef(k/2 + j_prime) + master_list[j_index-1][k].coef(j_prime - k/2)
                    # beta
                    
                    if j_prime <= (j_index - 1 - k)/2:
                        beta_j_prime = master_list[j_index-1][k-1].coef( (j_index + 1 - k)/2 + j_prime ) - master_list[j_index-1][k-1].coef( (j_index - 1 - k)/2 - j_prime )
                    elif j_prime == (j_index - k) / 2:
                        beta_j_prime = master_list[j_index-1][k-1].coef( j_index - k + 1/2 )
                    else:
                        beta_j_prime = master_list[j_index-1][k-1].coef( (j_index + 1 - k)/2 + j_prime ) + master_list[j_index-1][k-1].coef( j_prime - (j_index + 1 - k)/2 )
                    
                    final_sol = 0.5 * (alpha_j_prime + beta_j_prime)
                    if cur_coef != final_sol:
                        print(f"DISCREPANCY for j={cur_j}, k={k}, j'={j_prime}: found coef = {cur_coef}, recursive rel = {final_sol}")
            
            
        
        for k in range(1, int(j_index + 1)):
            
            max_j = max(master_list[j_index][k].nonzero_coefs.keys())
            
            cur_tab_row = []
            for j in np.arange(0, max_j+1/2, 1/2):
                cur_tab_row.append(master_list[j_index][k].coef(j))
            cur_max_length = len(cur_tab_row)
            if actual_max_length < cur_max_length:
                actual_max_length = cur_max_length
            
            decomposition_table.append(cur_tab_row)
            if int(2*cur_j) % 2 == 0:
                decomposition_table_row_names.append(str(int(cur_j)) + ", " + str(k))
            else:
                decomposition_table_row_names.append(str(int(2*cur_j)) + "/2, " + str(k))
        #print("    New max length =", actual_max_length)
        if cur_j != j_upper_cap:
            decomposition_table_row_names.append("-")
    
    # we trim everything
    for i in range(len(decomposition_table)):
        decomposition_table[i] = decomposition_table[i][:actual_max_length]
    
    print_table("c^j_k", decomposition_table_column_names(actual_max_length), decomposition_table_row_names, decomposition_table, omit_strings=["0"])
    

print_decomposition_table_from_recursion(11/2) #THIS IS SO FAST


