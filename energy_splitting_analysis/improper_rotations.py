

import numpy as np

# We only deal with proper and improper rotations here; no funky shit

inversion_matrix = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])

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

def angle_to_m(angle, max_q = 12):
    # nontrivial - a rational fraction approximator
    
    if np.round(angle, 10) == 0.0:
        return([0, 1])
    
    x = angle / (2.0 * np.pi)
    
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
    return(constrained_angle(reduce_frac([best_p, best_q])))
    

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


class ImproperRotation():
    
    # Static variables
    rounding_decimals = 10 # so that sin pi/2 = 1 etc
    
    # Each operation is defined as follows:
    #   1. a rotation about an AXIS (given by a unit vector) by a specified ANGLE (right-handed convention)
    #   2. If INVERSION == True, an inversion about [0, 0, 0]
    
    def __init__(self, axis, multiplicity, inversion):
        
        # multiplicity is a 2-element list [p, q], so that q is the foldedness of the rotation axis, and the rot angle is 2pi * p/q. Both are integer
        
        #print(np.linalg.norm(axis))
        self.axis = np.round(np.array(axis) / np.linalg.norm(axis), decimals = ImproperRotation.rounding_decimals)
        self.multiplicity = reduce_frac(multiplicity)
        self.inversion = inversion
    
    def cartesian_rep(self):
        # matrix as a list of rows
        res = np.zeros((3, 3))
        u_x = self.axis[0]
        u_y = self.axis[1]
        u_z = self.axis[2]
        t   = m_to_angle(self.multiplicity)
        print(t / (2.0 * np.pi))
        res[0][0] = np.cos(t) + u_x * u_x * (1.0 - np.cos(t))
        res[0][1] = u_x * u_y * (1.0 - np.cos(t)) - u_z * np.sin(t)
        res[0][2] = u_x * u_z * (1.0 - np.cos(t)) + u_y * np.sin(t)
        
        res[1][0] = u_y * u_x * (1.0 - np.cos(t)) + u_z * np.sin(t)
        res[1][1] = np.cos(t) + u_y * u_y * (1.0 - np.cos(t))
        res[1][2] = u_y * u_z * (1.0 - np.cos(t)) - u_x * np.sin(t)
        
        res[2][0] = u_z * u_x * (1.0 - np.cos(t)) - u_y * np.sin(t)
        res[2][1] = u_z * u_y * (1.0 - np.cos(t)) + u_x * np.sin(t)
        res[2][2] = np.cos(t) + u_z * u_z * (1.0 - np.cos(t))
        
        if self.inversion:
            res = np.matmul(res, inversion_matrix)
        
        res = np.round(res, decimals = ImproperRotation.rounding_decimals)
        
        return(np.matrix(res))
    
    def SU2_rep(self):
        # matrix as list of rows
        
        res  = np.zeros((2, 2), dtype = np.complex_)
        t    = m_to_angle(self.multiplicity)
        w    = np.cos(t / 2.0)
        norm = np.abs(np.sin(t / 2.0))
        u_x  = self.axis[0] * norm
        u_y  = self.axis[1] * norm
        u_z  = self.axis[2] * norm
        
        res[0][0] = complex(w, u_x)
        res[0][1] = complex(-u_y, u_z)
        res[1][0] = complex(u_y, u_z)
        res[1][1] = complex(w, -u_x)
        
        if self.inversion:
            # TODO FIGURE THIS OUT!!!!!!
            print("inversion ignored!!!")
        
        res = np.round(res, decimals = ImproperRotation.rounding_decimals)
        return(np.matrix(res))
        
    
    def __eq__(self, SO):
        
        if self.inversion != SO.inversion:
            return(False)
        
        # special cases: no rotation
        if self.multiplicity[0] == 0 and SO.multiplicity[0] == 0:
            return(True)
        
        if np.array_equal(self.axis, SO.axis) and ((self.multiplicity[0] == SO.multiplicity[0] and self.multiplicity[1] == SO.multiplicity[1]) or m_to_angle(self.multiplicity) == m_to_angle(SO.multiplicity)):
            return(True)
        
        if np.array_equal(self.axis, -SO.axis) and self.multiplicity[0] == inverse_angle(SO.multiplicity)[0] and self.multiplicity[1] == inverse_angle(SO.multiplicity)[1]:
            return(True) # inverting the axis == inverting the angle
        return(False)
        
        """if self.op_type == "rotation":
            if SO.op_type == "rotation":
                if self.axis == SO.axis and self.angle == SO.angle:
                    return(True)
            return(False)
        if self.op_type == "inversion":
            if SO.op_type == "inversion":
                if self.axis == SO.axis or self.axis == -SO.axis:
                    return(True)
            return(False)"""
    
    def __add__(self, SO):
        # P_R + P_S = P_R P_S
        
        # This method seems to have a problem for when the result is meant to be E - we add some catching exceptions
        
        # First - if axes are parallel, then we just add multiplicities
        if np.array_equal(self.axis, SO.axis):
            new_axis = self.axis
            new_multiplicity = add_multiplicities(self.multiplicity, SO.multiplicity)
        elif np.array_equal(self.axis, -SO.axis):
            new_axis = self.axis
            
            new_multiplicity = add_multiplicities(self.multiplicity, inverse_angle(SO.multiplicity))
        
        else:
            #https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula
            angle1 = m_to_angle(self.multiplicity)
            angle2 = m_to_angle(SO.multiplicity)
            
            a_1, b_1, c_1, d_1 = Euler_Rodriguez_encoder(self.axis, angle1)
            a_2, b_2, c_2, d_2 = Euler_Rodriguez_encoder(SO.axis, angle2)
            
            a = a_1 * a_2 - b_1 * b_2 - c_1 * c_2 - d_1 * d_2
            b = a_1 * b_2 + b_1 * a_2 - c_1 * d_2 + d_1 * c_2
            c = a_1 * c_2 + c_1 * a_2 - d_1 * b_2 + b_1 * d_2
            d = a_1 * d_2 + d_1 * a_2 - b_1 * c_2 + c_1 * b_2
            
            # here when a is supposed to equal 1 or -1, it often jumps over by a bit: we map these cases unto hard-coded ones and minus ones
            if a < -1.0 and a > -1.001:
                a = -1.0
            if a > 1.0 and a < 1.001:
                a = 1.0
            
            
            new_axis, new_angle = Euler_Rodriguez_decoder(a, b, c, d)
            
            new_multiplicity = angle_to_m(new_angle)
            #print("HA", new_angle, new_multiplicity)
            
        if (self.inversion == True and SO.inversion == True) or (self.inversion == False and SO.inversion == False):
            new_inversion = False
        else:
            new_inversion = True
        
        """
        # if new axis has zero norm, this has to be the identity matrix
        if np.round(np.linalg.norm(new_axis), ImproperRotation.rounding_decimals) == 0.0:
            if new_inversion == True:
                return(identity_inversion)
            if new_inversion == False:
                return(identity_rotation)"""
        
        # rounding here
        new_axis = np.round(new_axis, decimals = ImproperRotation.rounding_decimals)
        #new_angle = np.round(new_angle, decimals = ImproperRotation.rounding_decimals)
        return(ImproperRotation(new_axis, new_multiplicity, new_inversion))
    
    def inverse(self):
        return(ImproperRotation(self.axis, inverse_angle(self.multiplicity), self.inversion))
    
    def __neg__(self):
        return(self.inverse())
    
    def __sub__(self, SO):
        return(self + SO.inverse())
        

# ------------------- Common instances ---------------------

identity_rotation = ImproperRotation([1.0, 0.0, 0.0], [0, 1], False)
identity_inversion = ImproperRotation([1.0, 0.0, 0.0], [0, 1], True)


