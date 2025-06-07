import numpy as np
import matplotlib.pyplot as plt
import os
import statistics
import sys

data = np.array([99.0, 98.5, 97.5, 98.5, 97.5, 97.0, 98.0, 97.5, 97.0])
upstream_laser = -779.95
bac = np.array([156.5, 155.5, 156.5, 156.0, 156.0, 156.0])

htof = np.array([
    [343.5, 345.0, 343.5],  # U1
    [342.0, 344.0, 343.0],  # U2
    [341.0, 344.0, 343.0],  # U3
    [342.0, 345.0, 343.5],  # U4
])


htof2 = np.array([
    344.0, 342.0, 342.5, 344.0,
    346.0, 342.5, 343.0, 345.0,
    343.5, 341.5, 343.0, 343.0
])


sac = np.array([
    114.5, 112.5, 113.0,
    114.5, 112.5, 113.0,
    113.5, 112.5, 112.5
])

kvc2 = np.array([
    96.5, 95.5, 95.0,
    96.5, 96.0, 95.0,
    96.5, 96.0, 96.0
])

bh2 = np.array([
    224.0, 223.5, 222.5,
    222.0, 222.0, 221.0,
    224.0, 223.5, 223.0 
])

bh2_2 = np.array([
    402.0, 402.5, 402.0
])

htof_window_up = np.array([72.5, 73.5])
htof_window_down = np.array([37.0, 37.5])

print(upstream_laser + statistics.mean(kvc2) + 10.0)

sys.exit()

print(statistics.mean(htof_window_up))
print(statistics.mean(htof_window_down))
print(statistics.mean(data) + upstream_laser + 10.0)
print(-statistics.mean(bac) + upstream_laser)
print(-statistics.mean(bac) + upstream_laser + 1300)

print(-statistics.mean(htof.flatten()) - 150.071)
print(statistics.mean(bh2.flatten()) + upstream_laser)


# kvc1_data = {
#     # z, y
#     3: [28.875, 59.25],
#     4: [60.75,  60.75],
#     5: [61.5,   62.25],
#     6: [31.0,   57.75]
# }
kvc1_data = {
    # z, y
    3: [28.875+9.0, 59.25],
    4: [60.75+9.0,  60.75],
    5: [61.5+9.0,   62.25],
    6: [31.0+9.0,   57.75]
}



def theta1(up, down):
    return np.asin( (up[0]-down[0])/(up[1]+down[1]) ) * 180.0/np.pi

def z_kvc1cm(up, down):
    theta = np.asin( (up[0]-down[0])/(up[1]+down[1]) )
    z1 = (down[1]*(up[0]-down[0])) / ((up[1]+down[1])*np.cos(theta))
    z2 = down[0]/np.cos(theta)
    # return z1 + z2 + 696.141
    return z1 + z2

print(4, 3, theta1(kvc1_data[4], kvc1_data[3]))
print(4, 6, theta1(kvc1_data[4], kvc1_data[6]))

print(5, 3, theta1(kvc1_data[5], kvc1_data[3]))
print(5, 6, theta1(kvc1_data[5], kvc1_data[6]))


print(4, 3, z_kvc1cm(kvc1_data[4], kvc1_data[3]))
print(4, 6, z_kvc1cm(kvc1_data[4], kvc1_data[6]))

print(5, 3, z_kvc1cm(kvc1_data[5], kvc1_data[3]))
print(5, 6, z_kvc1cm(kvc1_data[5], kvc1_data[6]))

print(696.141 + statistics.mean([
    z_kvc1cm(kvc1_data[4], kvc1_data[3]),
    z_kvc1cm(kvc1_data[4], kvc1_data[6]),
    z_kvc1cm(kvc1_data[5], kvc1_data[3]),
    z_kvc1cm(kvc1_data[5], kvc1_data[6])
]
))

kvc1_data = {
    1: [243.5, 288.0, 0],
    2: [243.6, 289.3, 0],
    3: [243.5, 301.0, 53.2],
    4: [243.6, 302.5, 53.2],
    5: [243.5, 274.5, 53.5],
    6: [244.0, 276.2, 53.5],
}

# def theta2(up, down):
#     return np.atan( (up[1]-down[1])/(up[2]+down[2]) ) * 180.0/np.pi

# print(3, 5, theta2(kvc1_data[3], kvc1_data[5]))
# print(4, 6, theta2(kvc1_data[4], kvc1_data[6]))

# print(1, 5, theta2(kvc1_data[1], kvc1_data[5]))
# print(2, 6, theta2(kvc1_data[2], kvc1_data[6]))

# print(3, 1, theta2(kvc1_data[3], kvc1_data[1]))
# print(4, 2, theta2(kvc1_data[4], kvc1_data[2]))

# print(kvc1_data[1][1]-kvc1_data[1][0])
# print(kvc1_data[2][1]-kvc1_data[2][0])
