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

def rotate_list_lines(list_lines, point_on_axis, rot_axis, angle):
    list_rotated_lines = []
    for line in list_lines:
        cur_line_i = rotate_position(line[0], point_on_axis, rot_axis, angle)
        cur_line_f = rotate_position(line[1], point_on_axis, rot_axis, angle)
        list_rotated_lines.append([cur_line_i, cur_line_f])
    return(list_rotated_lines)

def project_position_onto_plane(position, point_on_plane, plane_normal_vector):
    unit_normal_vector = np.array(plane_normal_vector) / np.linalg.norm(plane_normal_vector)
    displacement = np.array(position) - np.array(point_on_plane)
    distance_from_plane = np.dot(displacement, unit_normal_vector)
    projected_position = np.array(position) - distance_from_plane * unit_normal_vector
    return(projected_position)

def planar_representation(position, new_basis):
    # position and point_on_plane are 3D vectors
    # basis is a list of 2 vectors
    new_basis_x = np.array(new_basis[0]) / np.linalg.norm(new_basis[0])
    new_basis_y = np.array(new_basis[1]) / np.linalg.norm(new_basis[1])
    planar_position = [np.dot(np.array(position), new_basis_x), np.dot(np.array(position), new_basis_y)]
    distance = np.dot(position, np.cross(new_basis_x, new_basis_y))
    
    return(planar_position, distance)

def list_planar_representation(list_position, new_basis, separate_by_distance = False):
    if separate_by_distance == False:
        list_planar_position = []
        for pos in list_position:
            cur_planar_position, cur_planar_distance = planar_representation(pos, new_basis)
            list_planar_position.append(cur_planar_position)
        return(list_planar_position)
    else:
        # return {dist: [list planar position]}
        plan_pos_by_dist = {}
        encountered_dist_vals = []
        for pos in list_position:
            cur_planar_position, cur_planar_distance = planar_representation(pos, new_basis)
            cur_planar_distance = np.round(cur_planar_distance, decimals = 5)
            if cur_planar_distance in encountered_dist_vals:
                plan_pos_by_dist[cur_planar_distance].append(cur_planar_position)
            else:
                plan_pos_by_dist[cur_planar_distance] = [cur_planar_position]
                encountered_dist_vals.append(cur_planar_distance)
        return(plan_pos_by_dist)
        

def list_lines_planar_representation(list_lines, new_basis):
    list_planar_lines = []
    for line in list_lines:
        cur_planar_pos_i, cur_planar_dist_i = planar_representation(line[0], new_basis)
        cur_planar_pos_f, cur_planar_dist_f = planar_representation(line[1], new_basis)
        list_planar_lines.append([cur_planar_pos_i, cur_planar_pos_f])
    return(list_planar_lines)

def create_unit_cell(a, relative_pos = np.array([0.0, 0.0, 0.0])):
    positions_a = []
    positions_b = []
    lines = [] # [[vec1, vec2], ...]
    
    # cube corners
    positions_a.append(np.array([0, 0, 0]) + relative_pos)
    positions_a.append(np.array([a, 0, 0]) + relative_pos)
    positions_a.append(np.array([0, a, 0]) + relative_pos)
    positions_a.append(np.array([0, 0, a]) + relative_pos)
    positions_a.append(np.array([a, a, 0]) + relative_pos)
    positions_a.append(np.array([0, a, a]) + relative_pos)
    positions_a.append(np.array([a, 0, a]) + relative_pos)
    positions_a.append(np.array([a, a, a]) + relative_pos)
    
    lines.append([np.array([0, 0, 0]) + relative_pos, np.array([a, 0, 0]) + relative_pos])
    lines.append([np.array([0, 0, 0]) + relative_pos, np.array([0, a, 0]) + relative_pos])
    lines.append([np.array([0, 0, 0]) + relative_pos, np.array([0, 0, a]) + relative_pos])
    lines.append([np.array([a, 0, 0]) + relative_pos, np.array([a, a, 0]) + relative_pos])
    lines.append([np.array([a, 0, 0]) + relative_pos, np.array([a, 0, a]) + relative_pos])
    lines.append([np.array([0, a, 0]) + relative_pos, np.array([a, a, 0]) + relative_pos])
    lines.append([np.array([0, a, 0]) + relative_pos, np.array([0, a, a]) + relative_pos])
    lines.append([np.array([0, 0, a]) + relative_pos, np.array([a, 0, a]) + relative_pos])
    lines.append([np.array([0, 0, a]) + relative_pos, np.array([0, a, a]) + relative_pos])
    lines.append([np.array([a, a, 0]) + relative_pos, np.array([a, a, a]) + relative_pos])
    lines.append([np.array([0, a, a]) + relative_pos, np.array([a, a, a]) + relative_pos])
    lines.append([np.array([a, 0, a]) + relative_pos, np.array([a, a, a]) + relative_pos])
    
    # face centres
    positions_a.append(np.array([a/2, a/2, 0]) + relative_pos)
    positions_a.append(np.array([a/2, a/2, a]) + relative_pos)
    positions_a.append(np.array([a/2, 0, a/2]) + relative_pos)
    positions_a.append(np.array([a/2, a, a/2]) + relative_pos)
    positions_a.append(np.array([0, a/2, a/2]) + relative_pos)
    positions_a.append(np.array([a, a/2, a/2]) + relative_pos)
    
    # second type atoms
    positions_b.append(np.array([a/4, a/4, a/4]) + relative_pos)
    positions_b.append(np.array([3 * a/4, 3 * a/4, a/4]) + relative_pos)
    positions_b.append(np.array([a/4, 3 * a/4, 3 * a/4]) + relative_pos)
    positions_b.append(np.array([3 * a/4, a/4, 3 * a/4]) + relative_pos)
    
    lines.append([np.array([0, 0, 0]) + relative_pos, np.array([a/4, a/4, a/4]) + relative_pos])
    lines.append([np.array([a/2, a/2, 0]) + relative_pos, np.array([a/4, a/4, a/4]) + relative_pos])
    lines.append([np.array([a/2, 0, a/2]) + relative_pos, np.array([a/4, a/4, a/4]) + relative_pos])
    lines.append([np.array([0, a/2, a/2]) + relative_pos, np.array([a/4, a/4, a/4]) + relative_pos])
    
    lines.append([np.array([a, 0, a]) + relative_pos, np.array([3 * a/4, a/4, 3 * a/4]) + relative_pos])
    lines.append([np.array([a/2, 0, a/2]) + relative_pos, np.array([3 * a/4, a/4, 3 * a/4]) + relative_pos])
    lines.append([np.array([a, a/2, a/2]) + relative_pos, np.array([3 * a/4, a/4, 3 * a/4]) + relative_pos])
    lines.append([np.array([a/2, a/2, a]) + relative_pos, np.array([3 * a/4, a/4, 3 * a/4]) + relative_pos])
    
    lines.append([np.array([0, a, a]) + relative_pos, np.array([a/4, 3 * a/4, 3 * a/4]) + relative_pos])
    lines.append([np.array([a/2, a, a/2]) + relative_pos, np.array([a/4, 3 * a/4, 3 * a/4]) + relative_pos])
    lines.append([np.array([0, a/2, a/2]) + relative_pos, np.array([a/4, 3 * a/4, 3 * a/4]) + relative_pos])
    lines.append([np.array([a/2, a/2, a]) + relative_pos, np.array([a/4, 3 * a/4, 3 * a/4]) + relative_pos])
    
    lines.append([np.array([a, a, 0]) + relative_pos, np.array([3 * a/4, 3 * a/4, a/4]) + relative_pos])
    lines.append([np.array([a/2, a/2, 0]) + relative_pos, np.array([3 * a/4, 3 * a/4, a/4]) + relative_pos])
    lines.append([np.array([a/2, a, a/2]) + relative_pos, np.array([3 * a/4, 3 * a/4, a/4]) + relative_pos])
    lines.append([np.array([a, a/2, a/2]) + relative_pos, np.array([3 * a/4, 3 * a/4, a/4]) + relative_pos])
    
    return(positions_a, positions_b, lines)

def create_reduced_unit_cell(a, relative_pos = np.array([0.0, 0.0, 0.0])):
    # this one is periodic without repeating faces and lines
    positions_a = []
    positions_b = []
    lines = [] # [[vec1, vec2], ...]
    
    # cube corners
    positions_a.append(np.array([0, 0, 0]) + relative_pos)
    
    lines.append([np.array([0, 0, 0]) + relative_pos, np.array([a, 0, 0]) + relative_pos])
    lines.append([np.array([0, 0, 0]) + relative_pos, np.array([0, a, 0]) + relative_pos])
    lines.append([np.array([0, 0, 0]) + relative_pos, np.array([0, 0, a]) + relative_pos])
    
    # face centres
    positions_a.append(np.array([a/2, a/2, 0]) + relative_pos)
    positions_a.append(np.array([a/2, 0, a/2]) + relative_pos)
    positions_a.append(np.array([0, a/2, a/2]) + relative_pos)
    
    # second type atoms
    positions_b.append(np.array([a/4, a/4, a/4]) + relative_pos)
    positions_b.append(np.array([3 * a/4, 3 * a/4, a/4]) + relative_pos)
    positions_b.append(np.array([a/4, 3 * a/4, 3 * a/4]) + relative_pos)
    positions_b.append(np.array([3 * a/4, a/4, 3 * a/4]) + relative_pos)
    
    lines.append([np.array([0, 0, 0]) + relative_pos, np.array([a/4, a/4, a/4]) + relative_pos])
    lines.append([np.array([a/2, a/2, 0]) + relative_pos, np.array([a/4, a/4, a/4]) + relative_pos])
    lines.append([np.array([a/2, 0, a/2]) + relative_pos, np.array([a/4, a/4, a/4]) + relative_pos])
    lines.append([np.array([0, a/2, a/2]) + relative_pos, np.array([a/4, a/4, a/4]) + relative_pos])
    
    lines.append([np.array([a, 0, a]) + relative_pos, np.array([3 * a/4, a/4, 3 * a/4]) + relative_pos])
    lines.append([np.array([a/2, 0, a/2]) + relative_pos, np.array([3 * a/4, a/4, 3 * a/4]) + relative_pos])
    lines.append([np.array([a, a/2, a/2]) + relative_pos, np.array([3 * a/4, a/4, 3 * a/4]) + relative_pos])
    lines.append([np.array([a/2, a/2, a]) + relative_pos, np.array([3 * a/4, a/4, 3 * a/4]) + relative_pos])
    
    lines.append([np.array([0, a, a]) + relative_pos, np.array([a/4, 3 * a/4, 3 * a/4]) + relative_pos])
    lines.append([np.array([a/2, a, a/2]) + relative_pos, np.array([a/4, 3 * a/4, 3 * a/4]) + relative_pos])
    lines.append([np.array([0, a/2, a/2]) + relative_pos, np.array([a/4, 3 * a/4, 3 * a/4]) + relative_pos])
    lines.append([np.array([a/2, a/2, a]) + relative_pos, np.array([a/4, 3 * a/4, 3 * a/4]) + relative_pos])
    
    lines.append([np.array([a, a, 0]) + relative_pos, np.array([3 * a/4, 3 * a/4, a/4]) + relative_pos])
    lines.append([np.array([a/2, a/2, 0]) + relative_pos, np.array([3 * a/4, 3 * a/4, a/4]) + relative_pos])
    lines.append([np.array([a/2, a, a/2]) + relative_pos, np.array([3 * a/4, 3 * a/4, a/4]) + relative_pos])
    lines.append([np.array([a, a/2, a/2]) + relative_pos, np.array([3 * a/4, 3 * a/4, a/4]) + relative_pos])
    
    return(positions_a, positions_b, lines)


def generate_FCC_lattice(a, N, corner_pos = False):
    # N unit cells in each direction
    total_a = []
    total_b = []
    total_lines = []
    
    if corner_pos == False:
        corner_pos = - a * N / 2 * np.array([1.0, 1.0, 1.0])
    
    delta_x = np.array([a, 0.0, 0.0])
    delta_y = np.array([0.0, a, 0.0])
    delta_z = np.array([0.0, 0.0, a])
    for i in range(N):
        for j in range(N):
            for k in range(N):
                cur_a, cur_b, cur_lines = create_reduced_unit_cell(a, corner_pos + delta_x * i + delta_y * j + delta_z * k)
                total_a += cur_a
                total_b += cur_b
                total_lines += cur_lines
    return(total_a, total_b, total_lines)


# ------------------ PLOTTING --------------------------
def plot_list_pos(list_pos, color = -1, size = 50):
    list_pos_x = []
    list_pos_y = []
    for pos in list_pos:
        list_pos_x.append(pos[0])
        list_pos_y.append(pos[1])
    if color == -1:
        plt.scatter(list_pos_x, list_pos_y, s = size)
    else:
        plt.scatter(list_pos_x, list_pos_y, color = color, s = size)

def plot_list_lines(list_lines, style = 'dashed', color = -1):
    for line in list_lines:
        if color == -1:
            plt.plot([line[0][0], line[1][0]], [line[0][1], line[1][1]], style)
        else:
            plt.plot([line[0][0], line[1][0]], [line[0][1], line[1][1]], style, color = color)




basis_001 = [[1, 0, 0], [0, 1, 0]]
basis_110 = [[-1, 1, 0], [0, 0, 1]]
#basis_111 = [[1, 1, 0], [-1, 1, 2]]
basis_111 = [[-1, 1, 0], [-1, -1, 2]] #THIS is 100% correct!!!!!

side_basis = [[1, 0.0, 0.2], [0.0, 1.0, 0.3]]

cur_basis = basis_111

"""
a, b, lines = create_unit_cell(1.0)#generate_FCC_lattice(1.0, 5)#create_unit_cell(1.0)
cur_rot_axis = np.array([1.0, 1.0, 1.0])
cur_rot_angle = np.pi / 3
a = rotate_list_position(a, np.array([0.0, 0.0, 0.0]), cur_rot_axis, cur_rot_angle)
b = rotate_list_position(b, np.array([0.0, 0.0, 0.0]), cur_rot_axis, cur_rot_angle)
lines = rotate_list_lines(lines, np.array([0.0, 0.0, 0.0]), cur_rot_axis, cur_rot_angle)
p_a = list_planar_representation(a, cur_basis)
p_b = list_planar_representation(b, cur_basis)
p_lines = list_lines_planar_representation(lines, cur_basis)

plot_list_pos(p_a, "blue")
plot_list_pos(p_b, "red")
plot_list_lines(p_lines, 'grey')
ax = plt.gca()
#ax.set_xlim(1, 3)
#ax.set_ylim(1, 3)
ax.set_aspect('equal', adjustable='box')
plt.show()"""


"""
a, b, lines = generate_FCC_lattice(1.0, 10)#create_unit_cell(1.0)
p_a = list_planar_representation(a, cur_basis, separate_by_distance = True)
p_b = list_planar_representation(b, cur_basis, separate_by_distance = True)

p_a_dist_vals = list(p_a.keys())
p_b_dist_vals = list(p_b.keys())
p_a_dist_vals.sort()
p_b_dist_vals.sort()

#print(p_a_dist_vals)
start_i = int(np.round(len(p_a_dist_vals)/2)-3)
for n in range(6):
    ax = plt.subplot(2, 3, n+1)
    plt.title("dist =" + str(p_a_dist_vals[start_i+n]))
    plot_list_pos(p_a[p_a_dist_vals[start_i+n]], "blue", size = 5)
    ax.set_xlim(4, 9)
    ax.set_ylim(2, 7)
    ax.set_aspect('equal', adjustable='box')

plt.tight_layout()
plt.show()

start_i = int(np.round(len(p_b_dist_vals)/2)-3)
for n in range(6):
    ax = plt.subplot(2, 3, n+1)
    plt.title("dist =" + str(p_b_dist_vals[start_i+n]))
    plot_list_pos(p_b[p_b_dist_vals[start_i+n]], "red", size = 5)
    ax.set_xlim(4, 9)
    ax.set_ylim(2, 7)
    ax.set_aspect('equal', adjustable='box')

plt.tight_layout()
plt.show()

# Conclusion: both planar section types are a) evenly spaced (sqrt(3) / 3), b) with in triplets (so moving 3 planes down returns back to the same position)
"""




a, b, lines = generate_FCC_lattice(1.0, 10)#create_unit_cell(1.0)
p_a = list_planar_representation(a, cur_basis, separate_by_distance = True)
p_b = list_planar_representation(b, cur_basis, separate_by_distance = True)

p_a_dist_vals = list(p_a.keys())
p_b_dist_vals = list(p_b.keys())
p_a_dist_vals.sort()
p_b_dist_vals.sort()

layer_colours = ["#cc3300", "#cc9900", "#99cc00"]

ax = plt.subplot(1, 2, 1)
start_i = int(np.round(len(p_a_dist_vals)/2)-1)
for n in range(3):
    plt.title("Normally oriented FCC")
    plt.plot([-2.0, 2.0], [0.0, 0.0], linestyle = "dashed", color = "grey")
    plt.plot([0.0, 0.0], [-2.0, 2.0], linestyle = "dashed", color = "grey")
    plot_list_pos(p_a[p_a_dist_vals[start_i+n]], layer_colours[n], size = 20)
    print("normal dist =", p_a_dist_vals[start_i+n])
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal', adjustable='box')

rot_axis_point = np.array([0.0, 0.0, 0.0])
rot_axis_vec = np.array([1.0, 1.0, 1.0])
rot_angle = 2.0 * np.pi / 6.0 #sixfold rotation

a_rot = rotate_list_position(a, rot_axis_point, rot_axis_vec, rot_angle)
b_rot = rotate_list_position(b, rot_axis_point, rot_axis_vec, rot_angle)

p_a_rot = list_planar_representation(a_rot, cur_basis, separate_by_distance = True)
p_b_rot = list_planar_representation(b_rot, cur_basis, separate_by_distance = True)

p_a_rot_dist_vals = list(p_a_rot.keys())
p_b_rot_dist_vals = list(p_b_rot.keys())
p_a_rot_dist_vals.sort()
p_b_rot_dist_vals.sort()

ax = plt.subplot(1, 2, 2)
start_i = int(np.round(len(p_a_rot_dist_vals)/2)-1)
for n in range(3):
    plt.title("FCC rotated by 1/6 of revol.")
    plt.plot([-2.0, 2.0], [0.0, 0.0], linestyle = "dashed", color = "grey")
    plt.plot([0.0, 0.0], [-2.0, 2.0], linestyle = "dashed", color = "grey")
    plot_list_pos(p_a_rot[p_a_rot_dist_vals[start_i+n]], layer_colours[n], size = 20)
    print("rotated dist =", p_a_rot_dist_vals[start_i+n])
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal', adjustable='box')



plt.tight_layout()
plt.show()

