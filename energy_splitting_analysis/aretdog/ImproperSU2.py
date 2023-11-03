

from .base_functions import *



class ImproperSU2():
    # an instance is an element (SU(2), inversion)
    # a * b = (SU(2)_a matmul SU_2(2)_b, i_a . i_b)
    
    def __init__(self, SU2_matrix, inversion):
        # SU2_matrix is an element of SU(2)
        # inversion is a boolean
        self.SU2_matrix = np.matrix(SU2_matrix)
        self.inversion = inversion
    
    def reverse(self):
        return(ImproperSU2(- self.SU2_matrix.copy(), self.inversion))
    
    
    def __add__(self, other):        
        return(ImproperSU2(np.matmul(self.SU2_matrix, other.SU2_matrix), self.inversion != other.inversion))
    def __mul__(self, other):
        return(self + other)
    
    def __eq__(self, other):
        return(np.all(self.SU2_matrix == other.SU2_matrix) and self.inversion == other.inversion)
    
    def find_closest_operation(self, improper_su2_operations):
        # returns the index of the operation closest to self
        i = -1
        smallest_difference = 1e9
        m_flat = np.array(self.SU2_matrix.flatten())
        for j in range(len(improper_su2_operations)):
            if improper_su2_operations[j].inversion != self.inversion:
                continue
            ar_dif = (m_flat - np.array(improper_su2_operations[j].SU2_matrix.flatten()))
            cur_dif = np.real(np.sum(ar_dif * np.conjugate(ar_dif)))
            if cur_dif < smallest_difference:
                smallest_difference = cur_dif
                i = j
        return(i)


