#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "definitions.h"
#include "io.h"
#include "kdtree.h"
#include "kdtools.h"

#define to_msec(a) a * 1000 / CLOCKS_PER_SEC

int main()
{
  // Couters
  int i, j, c, t, p, msec, mis = 0;
  clock_t pre_tree, dt_tree = 0, pre_tetra, dt_tetra, pre_mis, dt_mis;
  long total_time;

  // Initialisation
  //// Buffers allocation (should I put this in a struct??)
  int *idx, n;
  double  *bb1        = malloc(nsd*sizeof(double)),
          *bb2        = malloc(nsd*sizeof(double)),
          *point      = malloc(nsd*sizeof(double)),
         **mxyz_tetra = malloc(nsd*sizeof(double*)),
         **data_tetra = malloc(ndf*sizeof(double*)),
         **T          = malloc(nsd*sizeof(double*)),
         **TinvT      = malloc(nsd*sizeof(double*));
  for (i=0; i<nsd; i++){mxyz_tetra[i] = malloc(nen*sizeof(double));}
  for (i=0; i<ndf; i++){data_tetra[i] = malloc(nen*sizeof(double));}
  for (i=0; i<nsd; i++){T[i]          = malloc(nsd*sizeof(double));}
  for (i=0; i<nsd; i++){TinvT[i]      = malloc(nsd*sizeof(double));}
  KDTreeNode_t tree, leaf, *leaves = malloc(1000*sizeof(KDTreeNode_t));

  //// Read arrays
  double **mxyz_fine   = read_d2D(nsd, nn_fine,   "Data/mxyz_ns.fine");
  double **mxyz_coarse = read_d2D(nsd, nn_coarse, "Data/mxyz_ns.coarse");
  double **data_coarse = read_d2D(ndf, nn_coarse, "Data/data_ns.coarse");
  int    **mien_coarse = read_i2D(nen, ne_coarse, "Data/mien_ns.coarse");

  //// Output array
  double **data_fine = emptydArray(nsd, nn_fine);

  // Projection
  //// Parameters
  int    cap = 16;
  double tol = 0.01;

  //// Tree creation
  pre_tree = clock();
  tree = KDTree(mxyz_fine, cap, nn_fine);
  dt_tree += clock() - pre_tree;

  //// Tetrahedra Queries (Bounding Boxes)
  pre_tetra = clock();
  for (t=0; t<ne_coarse; t++) {
    nextTetra   (t, mien_coarse, mxyz_coarse, data_coarse, mxyz_tetra,
                 data_tetra, bb1, bb2, T, TinvT, tol);
    queryBBs    (tree, bb1, bb2, leaves, &c);
    projectTetra(leaves, &c, bb1, bb2, mxyz_tetra, data_tetra, data_fine, TinvT,
                 tol);
  }
  dt_tetra = clock() - pre_tetra;

  //// Missing Points Queries (Closest Neighbour)
  pre_mis = clock();
  for (p=0; p<nn_fine; p++) {
    if (isnan(data_fine[0][t])) {
      projectCN(tree, p, mxyz_fine, data_fine);
      mis += 1;
    }
  }
  dt_mis = clock() - pre_mis;
  total_time = to_msec(dt_tree) + to_msec(dt_tetra) + to_msec(dt_mis);

  // Printing
  printf("Cap           \t%d\n", cap);
  printf("Tree build    \t%ld milliseconds\n", to_msec(dt_tree));
  printf("BBox queries  \t%ld milliseconds\n", to_msec(dt_tetra));
  printf("%d CN queries \t%ld milliseconds\n", mis, to_msec(dt_mis));
  printf("TOTAL         \t%ld milliseconds\n", total_time);
  printf("\n");
  // printf("[%d]\t%ld\n", cap, to_msec(dt_tree)+to_msec(dt_tetra)+to_msec(dt_mis));

  return 0;
}
