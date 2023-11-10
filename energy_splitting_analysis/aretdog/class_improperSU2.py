

from .base_functions import *



class ImproperSU2():
    # an instance is an element (SU(2), inversion)
    # a * b = (SU(2)_a matmul SU_2(2)_b, i_a . i_b)
    
    def __init__(self, SU2_matrix, inversion, parent_improper_rotation):
        # SU2_matrix is an element of SU(2)
        # inversion is a boolean
        self.SU2_matrix = np.matrix(SU2_matrix)
        self.inversion = inversion
        self.parent_improper_rotation = parent_improper_rotation
    
    def reverse(self):
        return(ImproperSU2(- self.SU2_matrix.copy(), self.inversion, self.parent_improper_rotation))
    
    
    def __add__(self, other):        
        return(ImproperSU2(np.matmul(self.SU2_matrix, other.SU2_matrix), self.inversion != other.inversion, self.parent_improper_rotation + other.parent_improper_rotation))
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
    
    def wigner_d_matrix(self, j):
        dim = int(2*j+1)
        # First, we obtain the euler angles in the Z-Y-Z convention
        q_r = np.real(self.SU2_matrix.item(0,0))
        q_x = np.imag(self.SU2_matrix.item(0,0))
        q_y = np.real(self.SU2_matrix.item(1,0))
        q_z = np.imag(self.SU2_matrix.item(1,0))
        
        alpha = np.arctan2(q_z, q_r) - np.arctan2(-q_x, q_y)
        cos_beta = 2.0 * (q_r * q_r + q_z * q_z) - 1.0
        if cos_beta > 1.0:
            cos_beta = 1.0
        if cos_beta < -1.0:
            cos_beta = -1.0
        beta = np.arccos(cos_beta)
        gamma = np.arctan2(q_z, q_r) + np.arctan2(-q_x, q_y)
        
        # Using Wigner's formula for the small
        d_small = wigner_small_d_matrix(beta, j)
        D_matrix = np.zeros((dim, dim), dtype=complex)
        for a in range(dim):
            m_a = -j + a
            for b in range(dim):
                m_b = -j + b
                D_matrix[b][a] = np.exp(-1j * m_b * alpha) * d_small[b][a] * np.exp(-1j * m_a * gamma)
        # if j = 1, 3, 5 etc and this is inversion, flip the sign of the matrix. Just because. Dont ask.
        if int(j) == j:
            if j % 2 == 1 and self.inversion:
                D_matrix *= -1
        return(np.matrix(D_matrix))
    
    
    


