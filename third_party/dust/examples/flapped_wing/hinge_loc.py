import numpy as np

# Clear equivalent not needed in Python
sweep = 10
twist = 0
dihed = 5
n1 = np.array([0.5, 2, 0])
n2 = np.array([0.5, 4, 0])


def R(x, y, z):
    # Rotation matrices using numpy's deg2rad for degree to radian conversion
    Rz = np.array([[np.cos(np.deg2rad(z)), -np.sin(np.deg2rad(z)), 0],
                   [np.sin(np.deg2rad(z)), np.cos(np.deg2rad(z)), 0],
                   [0, 0, 1]])
    
    Ry = np.array([[np.cos(np.deg2rad(y)), 0, np.sin(np.deg2rad(y))],
                   [0, 1, 0],
                   [-np.sin(np.deg2rad(y)), 0, np.cos(np.deg2rad(y))]])
    
    Rx = np.array([[1, 0, 0],
                   [0, np.cos(np.deg2rad(x)), -np.sin(np.deg2rad(x))],
                   [0, np.sin(np.deg2rad(x)), np.cos(np.deg2rad(x))]])
    
    return Rz @ Ry @ Rx

rot = R(dihed, twist, -sweep)

Node1 = rot @ n1
Node2 = rot @ n2

print("Node1:", Node1)
print("Node2:", Node2)
print("n1:", n1)
print("n2:", n2)