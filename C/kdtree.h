#ifndef KDTREE_H_INCLUDED
#define KDTREE_H_INCLUDED

/* typedef struct xxxxx{} xxxx_t; */
typedef struct KDTreeNode{
  int     *indexes, // Array of indexes for the points
           depth,   // Node depth
           dim,     // Longest dimension
           n;       // Number of points in the leaf

  double **points,  // Coordinates in Fortran alignment
           pivot;   // Pivot dim coordinate.

  struct KDTreeNode *left, *rght; // Children
} KDTreeNode_t;

double **emptydArray(int d1, int d2);

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
);

void KDTreeSlidingBuild(
            KDTreeNode_t  *node,       // Node to be initialized
            double       **points,     // 3-by-n array (contiguous coords)
            int           *indexes,    // Array of indexes for the points
            int            depth,      // Node depth
            double        *bb1,        // Minimum bound point
            double        *bb2,        // Maximum bound point
            int            maxdepth,   // KDTree depth
            int            n           // Number of points
);

KDTreeNode_t KDTree(
            double **points,    // 3-by-n array (contiguous coords.)
            int      maxdepth,  // KDTree depth
            int      n
);

void projectCN(
            KDTreeNode_t  node, // Node to be queried
            int           idx,  // Point index
            double      **xyz,  // mxyz_fine
            double      **data  // data_fine
);

void queryBBs(
            KDTreeNode_t  node,       // Node to be queried<
            double       *bb1,        // Minimum point of the BB
            double       *bb2,        // Maximum point of the BB
            KDTreeNode_t *leaves_lst, // Out: List of found leaves
            int          *c           // Out: Found leaves count
);

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
);

#endif // KDTREE_H_INCLUDED
