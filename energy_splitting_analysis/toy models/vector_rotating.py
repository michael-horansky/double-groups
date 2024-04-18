import matplotlib.pyplot as plt
import numpy as np

def rotate_position(pos, point_on_axis, rot_axis, angle):
    a = np.array(pos) - np.array(point_on_axis)
    b = np.array(rot_axis) / np.linalg.norm(rot_axis)
    a_par_b = np.dot(a, b) * b
    a_perp_b = a - a_par_b
    if np.linalg.norm(a_perp_b) < 1e-05:
        # the vectors are parallel
        return(np.array(pos))
    omega = np.cross(b, a_perp_b)
    x_1 = np.cos(angle) / np.linalg.norm(a_perp_b)
    x_2 = np.sin(angle) / np.linalg.norm(omega)
    a_perp_b_rot = np.linalg.norm(a_perp_b) * (x_1 * a_perp_b + x_2 * omega)
    a_rot = a_perp_b_rot + a_par_b
    return(a_rot + np.array(point_on_axis))
    
def rotate_list_position(list_position, point_on_axis, rot_axis, angle):
    list_rotated_position = []
    for pos in list_position:
        list_rotated_position.append(rotate_position(pos, point_on_axis, rot_axis, angle))
    return(list_rotated_position)


y = np.array([0.0, 1.0, 0.0])
z = np.array([0.0, 0.0, 1.0])
m = np.array([1.0, -1.0, 0.0])
a = np.array([1.0, 1.0, 1.0])
rot_axis = np.array([1.0, -1.0, 0.0])
rot_angle = np.arccos(1/np.sqrt(3.0))

my_list = [y, z, m, a]
new_list = rotate_list_position(my_list, np.array([0.0, 0.0, 0.0]), rot_axis, rot_angle)
print("Cz_2 =", new_list[1])
print("Cy_2 =", new_list[0])
print("C_3 =", new_list[3])
print("m =", new_list[2])
