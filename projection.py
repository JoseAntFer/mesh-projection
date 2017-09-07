# -*- coding: utf-8 -*-
from time import perf_counter as time
from math import sqrt
import numpy as np
from scipy.spatial import cKDTree
from numba import autojit

def readBinaryMatrix(file_name, nrows, ncols, dtype=float, swapped=False):
    from struct import unpack as npck

    tp, size = ('d', 8) if dtype is float else ('i', 4)
    fmt = ('>' if swapped else '<') + tp
    array = []
    with open(file_name, 'rb') as fl:
        for i in range(nrows):
            array.append( [npck(fmt, fl.read(size)) for _ in range(ncols)] )
    return np.array( array )[:, :, 0]

@autojit
def nextTetra(xyz_coarse, pres_coarse, mien, tetra, press, T, TinvT):
    for i in range(4):
        tetra[i] = xyz_coarse[mien[i]]
        press[i] = pres_coarse[mien[i]]
    
    # Transformation matrix
    for i in range(3):
        T[i] = tetra[i+1] - tetra[0]
    
    # Longest radius for the ball that contains the tetra
    d1 = T[0,0]*T[0,0] + T[0,1]*T[0,1] + T[0,2]*T[0,2]
    d2 = T[1,0]*T[1,0] + T[1,1]*T[1,1] + T[1,2]*T[1,2]
    d3 = T[2,0]*T[2,0] + T[2,1]*T[2,1] + T[2,2]*T[2,2]
    radius = sqrt(max(d1, d2, d3))

    # Inverse Transpose
    T11_T22_minus_T21_T12 = T[1,1]*T[2,2] - T[2,1]*T[1,2]
    T01_T22_minus_T12_T20 = T[1,0]*T[2,2] - T[1,2]*T[2,0]
    T01_T21_minus_T11_T20 = T[1,0]*T[2,1] - T[1,1]*T[2,0]
    
    det = +T[0,0]*(T11_T22_minus_T21_T12) \
    		 -T[0,1]*(T01_T22_minus_T12_T20) \
    		 +T[0,2]*(T01_T21_minus_T11_T20)
    			  
    invdet = 1/det
    TinvT[0,0] =   T11_T22_minus_T21_T12*invdet
    TinvT[1,0] =  -(T[0,1]*T[2,2] - T[0,2]*T[2,1])*invdet
    TinvT[2,0] =   (T[0,1]*T[1,2] - T[0,2]*T[1,1])*invdet
    
    TinvT[0,1] =  -T01_T22_minus_T12_T20*invdet
    TinvT[1,1] =   (T[0,0]*T[2,2] - T[0,2]*T[2,0])*invdet
    TinvT[2,1] =  -(T[0,0]*T[1,2] - T[1,0]*T[0,2])*invdet
    
    TinvT[0,2] =   (T01_T21_minus_T11_T20)*invdet
    TinvT[1,2] =  -(T[0,0]*T[2,1] - T[2,0]*T[0,1])*invdet
    TinvT[2,2] =   (T[0,0]*T[1,1] - T[1,0]*T[0,1])*invdet
    
    return radius

@autojit
def projectIfInside(T, a, xyz_fine, ball_idx, press, pres_fine, qzero, qone, close):
    # Check with tolerance and project
    for idx in ball_idx:
        # Take the coordinates of the point
        xyz = xyz_fine[idx]

        # Transformation
        tras = xyz - a
        
        # Checks octants
        x = T[0,0]*tras[0] + T[0,1]*tras[1] + T[0,2]*tras[2]
        if (x > qzero):
            y = T[1,0]*tras[0] + T[1,1]*tras[1] + T[1,2]*tras[2]
            if (y > qzero):
                z = T[2,0]*tras[0] + T[2,1]*tras[1] + T[2,2]*tras[2]
                if (z > qzero):
                    # Check Plane
                    plane = 1 - x - y - z
                    if plane > qzero:
                        # Precompute squares
                        xx, yy, zz = x**2, y**2, z**2
                        
                        # Distance to (0, 0, 0)
                        da = sqrt( xx + yy + zz )
                        if da < close: 
                            pres_fine[idx] = press[0]
                            continue
                        da = 1/da
                        
                        # Distance to (1, 0, 0)
                        db = sqrt( (x-1)**2 + yy + zz )
                        if db < close: 
                            pres_fine[idx] = press[1]
                            continue
                        db = 1/db
        
                        # Distance to (0, 1, 0)
                        dc = sqrt( xx + (y-1)**2 + zz )
                        if dc < close: 
                            pres_fine[idx] = press[2]
                            continue
                        dc = 1/dc
        
                        # Distance to (0, 0, 1)
                        dd = sqrt( xx + yy + (z-1)**2 )
                        if dd < close: 
                            pres_fine[idx] = press[3]
                            continue
                        dd = 1/dd
                        
                        # Total
                        tot = da + db + dc + dd
                        
                        pres_fine[idx] = (press[0]*da + press[1]*db + \
                                          press[2]*dc + press[3]*dd)/tot
                                 
                                 
def project(mien_coarse, xyz_coarse, xyz_fine, pres_coarse, tol):
    t = time()

    # Precomputed parameters
    close = 0.01 # Percentage of closiness to consider the pressure is the same
    qzero = 0 - tol
    qone  = 1 + tol
    ne    = len(mien_coarse)
    
    # Allocate
    pres_fine = np.zeros(len(xyz_fine), dtype=float)*np.nan
    tetra = np.zeros((4, 3), dtype=float)
    press = np.zeros(4, dtype=float)
    T = np.zeros((3, 3), dtype=float)
    TinvT = np.zeros((3, 3), dtype=float)

    # Build a tree
    tree = cKDTree(xyz_fine)
    
    # Loop through the tetrahedra
    for m, mien in enumerate(mien_coarse):
        # Gather and preprocess data for next tetrahedron
        radius = nextTetra(xyz_coarse, pres_coarse, mien, tetra, press, T, TinvT)
        
        # Indices of the fine nodes in the ball
        ball_idx = tree.query_ball_point(tetra[0], radius+tol)
                
        # Check and project
        projectIfInside(TinvT, tetra[0], xyz_fine, ball_idx, press, pres_fine, qzero, qone, close)
        
        if m % int(ne/10) == 0:
            print('%.1f' % (m/ne*100) + '%')

    print('Project time: %.3f' % (time() - t))
    print('Points not classified: %d' % np.isnan(pres_fine).sum())
    return pres_fine
