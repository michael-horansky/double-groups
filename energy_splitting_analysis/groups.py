from improper_rotations import *



















x = ImproperRotation([0, 0, 1], np.pi, False)
y = ImproperRotation([0, 1, 0], np.pi, False)

print("x =", x.cartesian_rep())
print("y =", y.cartesian_rep())
z = x + y
print("x + y =", z.cartesian_rep())
print("x + y axis =", z.axis, "; x + y angle =", z.angle)

w = ImproperRotation([1, 0, 0], np.pi, False)

print(x == w)
print(z == w)
