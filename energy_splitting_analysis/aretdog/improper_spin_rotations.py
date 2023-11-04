

from .ImproperSU2 import *

# TODO figure this out!!!!!

class ImproperSpinRotation():
    # an instance is an element (axis, multiplicity, time reversal, inversion)
    # this is based on the fact we can uniquely construct time reversal so that it corresponds to quaternions with NEGATIVE a
    # So we can ascribe time reversality to an element of SU(2) by checking the sign of the real component of its (1, 1) element
    
    rounding_decimals = 10
    
    def __init__(self, axis, multiplicity = False, time_reversal = False, inversion = False):
        # SU2_matrix is an element of SU(2)
        # inversion is a boolean
        
        # can be initialized from an instance ImproperSU2
        
        if type(axis) == ImproperSU2 and multiplicity == False and time_reversal == False and inversion == False:
            
            a11 = axis.SU2_matrix.item((0, 0))
            a21 = axis.SU2_matrix.item((1, 0))
            
            w = np.round(np.real(a11))
            u_x = np.imag(a11)
            u_y = np.real(a21)
            u_z = np.imag(a21)
            new_axis = np.array([u_x, u_y, u_z])
            
            self.time_reversal = w < 0
            self.inversion = axis.inversion
            
            # if axis length == 0: sin(phi/2) == 0 - > phi == 0 or 2pi
            if np.round(u_x * u_x + u_y * u_y + u_z * u_z, decimals = ImproperSpinRotation.rounding_decimals) == 0.0:
                self.axis = [1.0, 0.0, 0.0]
                self.multiplicity = [0, 1]
            else:
                self.axis = new_axis
                # Our convention:
                #   if angle == pi: non-time reversed 
                
                t = 2.0 * np.arccos(w)
                new_multiplicity = angle_to_m(t)
                self.multiplicity = new_multiplicity
                
                """if w > 0:
                    # strictly positive - a non-time-reversed rotation. t lies on (0, pi)
                    """
        
        else:
            self.axis = np.round(np.array(axis) / np.linalg.norm(axis), decimals = ImproperSpinRotation.rounding_decimals)
            self.multiplicity = multiplicity
            self.time_reversal = time_reversal
            self.inversion = inversion
    
    def __str__(self):
        res = ""
        if self.multiplicity[0] == 0 or self.multiplicity[1] == 1:
            if self.inversion:
                res += "inversion"
            else:
                res += "identity rotation"
        else:
            res += f"{self.multiplicity[0]}-fold rotation as {seitz_notation(self.axis, self.multiplicity[1], self.inversion)}"
        if self.time_reversal:
            res += " with time reversal"
        return(res)
    
    
    def reverse(self):
        return(ImproperSpinRotation(self.axis, self.multiplicity, not self.time_reversal, self.inversion))
    
    def ImproperSU2_rep(self):
        # we need to constraint the angle into a (-pi, pi) interval, so that cos(a/2) is positive
        res  = np.zeros((2, 2), dtype = np.complex_)
        t    = m_to_angle(self.multiplicity)
        while t < - np.pi:
            t += 2.0 * np.pi
        while t >= np.pi:
            t -= 2.0 * np.pi
        
        w    = np.cos(t / 2.0)
        norm = np.abs(np.sin(t / 2.0)) #TODO check if the abs is necessary! check algebraically
        u_x  = self.axis[0] * norm
        u_y  = self.axis[1] * norm
        u_z  = self.axis[2] * norm
        
        res[0][0] = complex(w, u_x)
        res[0][1] = complex(-u_y, u_z)
        res[1][0] = complex(u_y, u_z)
        res[1][1] = complex(w, -u_x)
        
        if self.time_reversal:
            res = - res
        
        return(ImproperSU2(np.matrix(res), self.inversion))
    
    
    def __add__(self, other):
        s1 = self.ImproperSU2_rep()
        s2 = other.ImproperSU2_rep()
        return(ImproperSpinRotation(s1 + s2))
    def __mul__(self, other):
        return(self + other)
    
    def __eq__(self, other):
        if self.time_reversal != other.time_reversal:
            return(False)
        if self.inversion != other.inversion:
            return(False)
        
        #same as for ImproperRotation
        if np.array_equal(self.axis, other.axis) and ((self.multiplicity[0] == other.multiplicity[0] and self.multiplicity[1] == other.multiplicity[1]) or m_to_angle(self.multiplicity) == m_to_angle(other.multiplicity)):
            return(True)
        
        #if np.array_equal(self.axis, -other.axis) and self.multiplicity[0] == inverse_angle(other.multiplicity)[0] and self.multiplicity[1] == inverse_angle(other.multiplicity)[1]:
        #    return(True) # inverting the axis == inverting the angle
        return(False)
    
        #return(self.ImproperSU2_rep() == other.ImproperSU2_rep())
    
    def find_closest_operation(self, improper_spin_rotations):
        # returns the index of the operation equal to self
        for i in range(len(improper_spin_rotations)):
            if improper_spin_rotations[i] == self:
                return(i)
        return(-1)


