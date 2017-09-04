# MESH PROJECTION ALGORITHM
#
# This version uses the KDTree data structure provided by scipy.spatial that
# uses balls to query enclosing subregions.
#
# The code is accelerated with Numba just in time compiler.
#
# This code is entirely sequential.
#
# Only the projection of one degree of freedom supported in this version.
#
# Notes: This code is not generating the exact same result than the one 
# provided in pres.fine although it is close. For tol around 1.2 it gets the 
# same number of points projected, so it is actually catching all the points as
# the other code does.

from projection import readBinaryMatrix, project

# Parameters
tol = 0.0
nn = 25014 ; nn2 = 110618
ne = 123911; ne2 = 570520

# Load data
print('Loading data...')
mien_coarse = readBinaryMatrix('mien.coarse', ne, 4, int, swapped=True)
xyz_coarse = readBinaryMatrix('mxyz.coarse', nn, 3, float, swapped=True)
xyz_fine = readBinaryMatrix('mxyz.fine', nn2, 3, float, swapped=True)
pres_coarse = readBinaryMatrix('pres.coarse', nn, 1, swapped=True).flatten()
mien_coarse -= 1 # Index in python starts in 0
print('Data loaded!')

# Projection
print('Projecting...')
pres_fine = project(mien_coarse, xyz_coarse, xyz_fine, pres_coarse, tol)

