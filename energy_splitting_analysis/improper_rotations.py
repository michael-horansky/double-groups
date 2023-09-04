

import numpy as np

# We only deal with proper and improper rotations here; no funky shit

inversion_matrix = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])

def Euler_Rodriguez_encoder(axis, angle):
    a = np.cos(angle / 2.0)
    b = axis[0] * np.sin(angle / 2.0)
    c = axis[1] * np.sin(angle / 2.0)
    d = axis[2] * np.sin(angle / 2.0)
    return(a, b, c, d)

def Euler_Rodriguez_decoder(a, b, c, d):
    axis = np.zeros(3)
    angle = 2.0 * np.arccos(a)
    sin_angle_half = np.sqrt(1.0 - a * a)
    axis[0] = b / sin_angle_half
    axis[1] = c / sin_angle_half
    axis[2] = d / sin_angle_half
    return(axis, angle)


class ImproperRotation():
    
    # Each operation is defined as follows:
    #   1. a rotation about an AXIS (given by a unit vector) by a specified ANGLE (right-handed convention)
    #   2. If INVERSION == True, an inversion about [0, 0, 0]
    
    def __init__(self, axis, angle, inversion):
        
        self.axis = np.array(axis) / np.linalg.norm(axis)
        self.angle = angle
        self.inversion = inversion
    
    def cartesian_rep(self):
        # matrix as a list of rows
        res = np.zeros((3, 3))
        u_x = self.axis[0]
        u_y = self.axis[1]
        u_z = self.axis[2]
        t   = self.angle
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
        
        """if self.op_type == "inversion":
            #https://en.wikipedia.org/wiki/Transformation_matrix#Reflection_2
            u_x = self.axis[0]
            u_y = self.axis[1]
            u_z = self.axis[2]
            res[0][0] = 1.0 - 2.0 * u_x * u_x
            res[0][1] = - 2.0 * u_x * u_y
            res[0][2] = - 2.0 * u_x * u_z
            
            res[1][0] = - 2.0 * u_x * u_y
            res[1][1] = 1.0 - 2.0 * u_y * u_y
            res[1][2] = - 2.0 * u_y * u_z
            
            res[2][0] = - 2.0 * u_x * u_z
            res[2][1] = - 2.0 * u_y * u_z
            res[2][2] = 1.0 - 2.0 * u_z * u_z"""
        
        res = np.round(res, decimals = 10) # so that sin 90 = 1 etc
        
        return(res)
            
        
    
    def __eq__(self, SO):
        
        if self.inversion != SO.inversion:
            return(False)
        if np.array_equal(self.axis, SO.axis) and self.angle == SO.angle:
            return(True)
        #if self.axis == -SO.axis and self.angle == -SO.angle:
        #    return(True)
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
        
        #https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula
        a_1, b_1, c_1, d_1 = Euler_Rodriguez_encoder(self.axis, self.angle)
        a_2, b_2, c_2, d_2 = Euler_Rodriguez_encoder(SO.axis, SO.angle)
        
        a = a_1 * a_2 - b_1 * b_2 - c_1 * c_2 - d_1 * d_2
        b = a_1 * b_2 + b_1 * a_2 - c_1 * d_2 + d_1 * c_2
        c = a_1 * c_2 + c_1 * a_2 - d_1 * b_2 + b_1 * d_2
        d = a_1 * d_2 + d_1 * a_2 - b_1 * c_2 + c_1 * b_2
        
        new_axis, new_angle = Euler_Rodriguez_decoder(a, b, c, d)
            
        if (self.inversion == True and SO.inversion == True) or (self.inversion == False and SO.inversion == False):
            new_inversion = False
        else:
            new_inversion = True
        
        # rounding here
        new_axis = np.round(new_axis, decimals = 10)
        return(ImproperRotation(new_axis, new_angle, new_inversion))
        




