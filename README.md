# affine
This is a script for calculating the affine coefficients (taking acount of distortion inaccuracy).

Two cameras are capturing the images of night sky with stars. 
Since the stars are far away (almost infinity, hehe), we may suppose that the scene is the same, and the only difference is 
in linear displacements, rotation (and distortion inaccuracy).

The user chooses the pairs of corresponding stars, get their coordinates (x1, y1) , (x2, y2).
Several points are needed (at least 4), because we need to solve a system of linear equations 
with the 7 unknown parameters = 6 (affine coefficients) + 1 (distortion correction coeff. of second[right] camera ).

The solution is approximate.
