"""
This is a script for calculating the affine coefficients (taking acount of distortion inaccuracy).

Two cameras are capturing the images of night sky with stars. 
Since the stars are far away (almost infinity, hehe), we may suppose that the scenes are the same, and the only difference is 
in linear displacements, rotation (and distortion inaccuracy).

The user chooses the pairs of corresponding stars, get their coordinates (x1, y1) , (x2, y2).
Several points are needed (at least 4 pairs), because we need to solve a system of linear equations 
with the 7 unknown parameters = 6 (affine coefficients) + 1 (distortion correction coeff. of second[right] camera ).
"""

import sys
from PIL import Image, ImageDraw
import numpy as np
np.set_printoptions(precision=2, suppress=True)

# Input of star pairs coordinates: (x, y) [at least 4 pairs]
width, height = 640, 480 # size of a picture 

fname = 'data/coords_640x480.txt'
coords = np.loadtxt(fname)

lX = coords[:, 0]
lY = coords[:, 1]
rX = coords[:, 2]
rY = coords[:, 3]

N = coords.shape[0] # number of pairs of points
M = coords.shape[1] # lX, lY, rX, rY == 4

if N < 4:
    print 'N < 4. Put some more pairs of points to file!'
    sys.exit(-1)

print 'Input coordinates from %r: \n' % fname, coords, '\n'



# Shift the coordinates: x = x - x_c; y = y - y_c
# So, now (0,0) point is in the center of image
# It is needed for distortion correction coefficient
x_c, y_c = width // 2, height // 2 # center coordinates
print 'x_center, y_center:', x_c, y_c, '\n'

lX -= x_c
rX -= x_c
lY -= y_c
rY -= y_c

print 'Shifted coordinates (with (0,0) in the center of image):'
print 'lX:', lX
print 'lY:', lY, '\n'
print 'rX', rX
print 'rY', rY, '\n'



# Solve the system of linear equations via pseudo inversion

# Below lX, lY, rX, rY mean lX[i], lY[i], rX[i], rY[i]:
# rX = a*lX + b*lY + e + e_2*rX*(rX**2 + rY**2)**0.5
# rY = c*lX + d*lY + f + e_2*rY*(rX**2 + rY**2)**0.5

# xi = A*f + nu,      where:

# xi.T = rX[0], rY[0] ... rX[N-1], rY[N-1]
# f.T = a, b, c, d, e, f, e_2  -- 6 affine coeff-s + distortion coeff
# nu -- vector of inaccuracy (?) I don't use it right now

# A = [ # L = N-1   -- Last 
# lX[0] lY[0]  0     0      1  0  rX[0]*(rX[0]**2 + rY[0]**2)**0.5 ;
# 0     0      lX[0] lY[0]  0  1  rY[0]*(rX[0]**2 + rY[0]**2)**0.5 ;
# ...                                                          ... ;
# ...                                                          ... ;
# ...                                                          ... ;
# lX[L] lY[L]  0     0      1  0  rX[L]*(rX[L]**2 + rY[L]**2)**0.5 ;
# 0     0      lX[L] lY[L]  0  1  rY[L]*(rX[L]**2 + rY[L]**2)**0.5 . ]

xi = np.zeros((2*N, 1))
f = np.zeros((7, 1))
a = np.zeros((2*N, 7)) # matrix A

for i in xrange(N): # fill the xi vector
    xi[2*i] = rX[i]
    xi[2*i + 1] = rY[i]

for i in xrange(N): # fill the A matrix
    tmp = rX[i]*(rX[i]**2 + rY[i]**2)**0.5
    a[2*i] = [lX[i], lY[i], 0, 0, 1, 0, tmp]
    
    tmp = rY[i]*(rX[i]**2 + rY[i]**2)**0.5
    a[2*i + 1] = [0, 0, lX[i], lY[i], 0, 1, tmp]

print 'xi = A*f + nu'
print 'xi.T:\n', xi.T
print 'A:\n', a, '\n'


pinv_a = np.linalg.pinv(a)

print 'Pseudo-inverted A -- pinv_a:\n', pinv_a, '\n'

f = np.dot(pinv_a, xi)

# Print the results
np.set_printoptions(precision=8, suppress=True)
print 'Result:'
print 'f.T:', f.T

# Test calculated coefficients applying them to clouds/stars images



