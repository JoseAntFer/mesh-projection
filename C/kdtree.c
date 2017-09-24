#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "kdtools.h"
#include "kdtree.h"
#include "definitions.h"

double **emptydArray(int d1, int d2)
/* Allocates an array full of NAN values */
{
  int i, j;
  double **A = malloc(d1*sizeof(double*));
  for (i=0; i<d1; i++) { A[i] = malloc(d2*sizeof(double)); }

  for (i=0; i<d1; i++) {
    for (j=0; j<d2; j++) {
      A[i][j] = NAN;
    }
  }

  return A;
}

void KDTreeSlidingBuild(
            KDTreeNode_t  *node,       // Node to be initialized
            double       **points,     // 3-by-n array (contiguous coords)
            int           *indexes,    // Array of indexes for the points
            int            depth,      // Node depth
            double        *bb1,        // Minimum bound point
            double        *bb2,        // Maximum bound point
            int            cap,        // KDTree depth
            int            n           // Number of points
                      )
{

    // Common to nodes and leaves
    node->depth = depth;
    node->dim   = longestDimension(bb1, bb2);

    // I am a node
    if (n > cap) {

        // I don't have points
        node->n = -1;

        // Compute partition
        node->pivot = 0.5*(bb1[node->dim] + bb2[node->dim]);

        //// Allocate data for the children
        int     *left_idx = malloc(sizeof(int)*n);
        int     *rght_idx = malloc(sizeof(int)*n);
        double **left_pts = malloc(sizeof(double*)*3);
        double **rght_pts = malloc(sizeof(double*)*3);
        left_pts[0] = malloc(sizeof(double)*n);
        left_pts[1] = malloc(sizeof(double)*n);
        left_pts[2] = malloc(sizeof(double)*n);
        rght_pts[0] = malloc(sizeof(double)*n);
        rght_pts[1] = malloc(sizeof(double)*n);
        rght_pts[2] = malloc(sizeof(double)*n);

        //// Classify with Slide
        int length[2] = {0, 0};
        node->pivot = classifyWithSlide(points, indexes, left_idx, rght_idx,
                          left_pts, rght_pts, node->pivot, node->dim, length, n);

        //// Children boundaries
        double left_bb2[3] = {bb2[0], bb2[1], bb2[2]},
               rght_bb1[3] = {bb1[0], bb1[1], bb1[2]};
        left_bb2[node->dim] = node->pivot;
        rght_bb1[node->dim] = node->pivot;

        //// Children creation
        ////// Left Child
        KDTreeNode_t *left_child = malloc(sizeof(KDTreeNode_t));
        KDTreeSlidingBuild(left_child, left_pts, left_idx, depth + 1,
                           bb1, left_bb2, cap, length[0]);
        node->left = left_child;

        ////// Right Child
        KDTreeNode_t *rght_child = malloc(sizeof(KDTreeNode_t));
        KDTreeSlidingBuild(rght_child, rght_pts, rght_idx, depth + 1,
                           rght_bb1, bb2, cap, length[1]);
        node->rght = rght_child;

    // I am a leaf
    } else {
        node->n       = n;
        node->indexes = indexes;
        node->points  = points;
    }
}

KDTreeNode_t KDTree(
                    double **points,    // 3-by-n array (contiguous coords.)
                    int      maxdepth,  // KDTree depth
                    int      n          // Length of the points given.
                   )
{
    int i;
    // Indexes vector
    int *indexes = malloc(sizeof(int)*n);
    for (i=0; i<nn_fine; i++) { indexes[i] = i; }

    // Bounding Box of the points
    double bb1[3], bb2[3];
    boundingBox(points, bb1, bb2, n, 0.0);

    // Initialization of the root
    KDTreeNode_t *root = malloc(sizeof(KDTreeNode_t));
    KDTreeSlidingBuild(root, points, indexes, 0, bb1, bb2, maxdepth, n);

    return *root;
}

void nextTetra(
               int      t,           // Tetra index to gather
               int    **mien_coarse, // Coarse connectivity array
               double **mxyz_coarse, // Coarse coordinates array
               double **data_coarse, // Coarse data array
               double **mxyz_tetra,  // Buffer tetrahedron coords array
               double **data_tetra,  // Buffer tetrahedron data array
               double  *bb1,         // Minimum bound point
               double  *bb2,         // Maximum bound point
               double **T,           // Buffer for Translation
               double **TinvT,       // Buffer for inv(Translation).T
               double   tol          // Tolerance for the BB
             )
{
/* Gather and preprocess information for the m tetrahedron in mien_coarse.
T and TinvT are C contiguous. */
  int i, j;

  // Gather tetra_xyz (mien_coarse can be transposed for clearliness!!)
  for (i=0; i<nen; i++) {
    for (j=0; j<nsd; j++) { mxyz_tetra[j][i] = mxyz_coarse[j][mien_coarse[i][t]-1]; }
    for (j=0; j<ndf; j++) { data_tetra[j][i] = data_coarse[j][mien_coarse[i][t]-1]; }
  }

  // Bounding box
  boundingBox(mxyz_tetra, bb1, bb2, nen, tol);

  // Transformation matrix (C contiguous)
  for (i=0; i<nsd; i++) {
    for (j=0; j<nsd; j++) {
      T[j][i] = mxyz_tetra[i][j+1] - mxyz_tetra[i][0];
    }
  }

  // Inverse Transpose
  double T11_T22_minus_T21_T12 = T[1][1]*T[2][2] - T[2][1]*T[1][2];
  double T01_T22_minus_T12_T20 = T[1][0]*T[2][2] - T[1][2]*T[2][0];
  double T01_T21_minus_T11_T20 = T[1][0]*T[2][1] - T[1][1]*T[2][0];

  double invdet = 1/( +T[0][0]*(T11_T22_minus_T21_T12) \
                      -T[0][1]*(T01_T22_minus_T12_T20) \
                      +T[0][2]*(T01_T21_minus_T11_T20) );

  TinvT[0][0] =   T11_T22_minus_T21_T12*invdet;
  TinvT[1][0] =  -(T[0][1]*T[2][2] - T[0][2]*T[2][1])*invdet;
  TinvT[2][0] =   (T[0][1]*T[1][2] - T[0][2]*T[1][1])*invdet;

  TinvT[0][1] =  -T01_T22_minus_T12_T20*invdet;
  TinvT[1][1] =   (T[0][0]*T[2][2] - T[0][2]*T[2][0])*invdet;
  TinvT[2][1] =  -(T[0][0]*T[1][2] - T[1][0]*T[0][2])*invdet;

  TinvT[0][2] =   (T01_T21_minus_T11_T20)*invdet;
  TinvT[1][2] =  -(T[0][0]*T[2][1] - T[2][0]*T[0][1])*invdet;
  TinvT[2][2] =   (T[0][0]*T[1][1] - T[1][0]*T[0][1])*invdet;
}

void queryBBs(
               KDTreeNode_t  node,       // Node to be queried<
               double       *bb1,        // Minimum point of the BB
               double       *bb2,        // Maximum point of the BB
               KDTreeNode_t *leaves_lst, // Out: List of found leaves
               int          *c           // Out: Found leaves count
              )
{
  // I am a leaf
  if (node.n >= 0) {
    leaves_lst[*c] = node;
    *c += 1;

  // I am a node
  } else {
    if (bb1[node.dim] < node.pivot) {
      if (bb2[node.dim] < node.pivot) {
        // Coincide both in left
        queryBBs(*node.left, bb1, bb2, leaves_lst, c);
      } else {
        // bb1 in left and bb2 in right
        double left_bb2[3] = {bb2[0], bb2[1], bb2[2]},
               rght_bb1[3] = {bb1[0], bb1[1], bb1[2]};
        left_bb2[node.dim] = node.pivot - 1e-10;
        rght_bb1[node.dim] = node.pivot;

        queryBBs(*node.left, bb1,      left_bb2, leaves_lst, c);
        queryBBs(*node.rght, rght_bb1, bb2,      leaves_lst, c);
      }
    } else {
      // Coincide both in right
      queryBBs(*node.rght, bb1, bb2, leaves_lst, c);
    }
  }
}

void projectTetra(
                  KDTreeNode_t *leaves,     // List of leaves to project with
                  int          *c,          // Number of leaves
                  double       *bb1,        // Minimum point of the BB
                  double       *bb2,        // Maximum point of the BB
                  double      **mxyz_tetra, // xyz coords of the tetrahadron
                  double      **data_tetra, // DoF for the tetrahedron nodes
                  double      **data_fine,  // DoF for the fine nodes
                  double      **TinvT,      // inv(Translation).T
                  double        tol         // Tolerance to project with
                )
{
  int l, p, d, idx;
  double px, py, pz, tras[3], xhi, eta, zeta, plane, xx, yy, zz, da, db, dc, dd, tot;
  double qone  = 1 + tol;
  double qzero = 0 - tol;
  double close = 0.05;

  // Leaves loop
  for (l=0; l<*c; l++) {
    // Points loop
    for (p=0; p<leaves[l].n; p++) {

      // Gather point data
      idx = leaves[l].indexes[p];
      /* I CAN PUT DATA_FINE CONTIGUOUSLY IN THE NODES AS WELL AS XYZ */
      if (!isnan(data_fine[0][idx])) { continue; }

      px = leaves[l].points[0][p];
      py = leaves[l].points[1][p];
      pz = leaves[l].points[2][p];

      // Check if inside
      //// Bounding box
      if ((px < bb1[0] - tol) || (px > bb2[0] + tol)) { continue; }
      if ((py < bb1[1] - tol) || (py > bb2[1] + tol)) { continue; }
      if ((pz < bb1[2] - tol) || (pz > bb2[2] + tol)) { continue; }

      //// Translation
      tras[0] = px - mxyz_tetra[0][0];
      tras[1] = py - mxyz_tetra[1][0];
      tras[2] = pz - mxyz_tetra[2][0];

      //// Checks octants
      xhi  = TinvT[0][0]*tras[0] + TinvT[0][1]*tras[1] + TinvT[0][2]*tras[2];
      if (xhi <= qzero) { continue; }
      eta  = TinvT[1][0]*tras[0] + TinvT[1][1]*tras[1] + TinvT[1][2]*tras[2];
      if (eta <= qzero) { continue; }
      zeta = TinvT[2][0]*tras[0] + TinvT[2][1]*tras[1] + TinvT[2][2]*tras[2];
      if (zeta <= qzero) { continue; }

      //// Check Plane
      plane = 1 - xhi - eta - zeta;
      if (plane <= qzero) { continue; }

      // Projection
      //// Precompute squares
      xx = px*px;   yy = py*py;   zz = pz*pz;

      //// Distance to (0, 0, 0)
      da = sqrt( xx + yy + zz );
      if (da < close) {
        for (d=0; d<ndf; d++){ data_fine[d][idx] = data_tetra[d][0]; }
        continue;
      }
      da = 1/da;

      //// Distance to (1, 0, 0)
      db = sqrt( (xhi-1)*(xhi-1) + yy + zz );
      if (db < close) {
        for (d=0; d<ndf; d++){ data_fine[d][idx] = data_tetra[d][1]; }
        continue;
      }
      db = 1/db;

      //// Distance to (0, 1, 0)
      dc = sqrt( xx + (eta-1)*(eta-1) + zz );
      if (dc < close) {
        for (d=0; d<ndf; d++){ data_fine[d][idx] = data_tetra[d][2]; }
        continue;
      }
      dc = 1/dc;

      //// Distance to (0, 0, 1)
      dd = sqrt( xx + yy + (zeta-1)*(zeta-1) );
      if (dd < close) {
        data_fine[idx] = data_tetra[3];
        continue;
      }
      dd = 1/dd;

      //// Total
      tot = 1/(da + db + dc + dd);
      for (d=0; d<ndf; d++){
        data_fine[d][idx] = (data_tetra[d][0]*da + \
                             data_tetra[d][1]*db + \
                             data_tetra[d][2]*dc + \
                             data_tetra[d][3]*dd)*tot;
      }
    }
  }

  *c = 0;
}

void projectCN(
          KDTreeNode_t  node, // Node to be queried
          int           idx,  // Point index
          double      **xyz,  // mxyz_fine
          double      **data  // data_fine
         )
{
  // I am a leaf
  if (node.n >= 0) {
    int p, best, global;
    double min_dist = INFINITY, dist, a, b, c;

    // Loop through the points in the leaf
    for (p=0; p<node.n; p++) {

      // No valid Neighbour (NAN)
      global = node.indexes[p];
      if (isnan(data[0][global])) { continue; } // Indexes into the struct

      // Squared distance computation
      a = node.points[0][p] - xyz[0][idx];
      b = node.points[1][p] - xyz[1][idx];
      c = node.points[2][p] - xyz[2][idx];
      dist = a*a + b*b + c*c;

      if (dist < min_dist) { best = global; }
    }

    // Copy the data
    for (p=0; p<nsd; p++) {
      data[p][idx] = data[p][best];
    }

  // I am a node
  } else {
    if (xyz[node.dim][idx] < node.pivot) {projectCN(*node.left, idx, xyz, data);}
    else                                 {projectCN(*node.rght, idx, xyz, data);}
  }
}
