#include <stdlib.h>
#include <stdio.h> // debugging
#include <math.h>

double maxd(double *array, int n)
/* Maximum value in a double array */
{
  int i;
  double maximum = array[0];

  for (i=1; i<n; i++) {
    if (array[i] > maximum) {
       maximum  = array[i];
    }
  }

  return maximum;
}


double mind(double *array, int n)
/* Maximum value in a double array */
{
  int i;
  double minimum = array[0];

  for (i=1; i<n; i++) {
    if (array[i] < minimum) {
       minimum  = array[i];
    }
  }

  return minimum;
}

void boundingBox(double **points, double *bb1, double *bb2, int n, double tol)
/* Stores the minimum bound point of array
in bb1 and the maximum in bb2 */
{
    bb1[0] = mind(points[0], n) - tol;
    bb1[1] = mind(points[1], n) - tol;
    bb1[2] = mind(points[2], n) - tol;

    bb2[0] = maxd(points[0], n) + tol;
    bb2[1] = maxd(points[1], n) + tol;
    bb2[2] = maxd(points[2], n) + tol;
}

int longestDimension(double *bb1, double *bb2)
/* Returns the index of the longest dimension */
{
    double best = fabs(bb2[0] - bb1[0]), current = 0.0;
    int    dim = 0, i;

    for (i=1; i<3; i++) {
      current = fabs(bb2[i] - bb1[i]);
        if (current > best) {
            dim = i;
            best = current;
        }
    }

    return dim;
}

void classify(double **points, int *indexes,
              int *left_idx, int *rght_idx,
              double **left_pts, double **rght_pts,
              double pivot, int dim, int *length, int n)
/* Fill the left and right points and indexes arrays using the pivot*/
{
    int lc = 0, rc = 0, i;

    for (i=0; i<n; i++) {
        if (points[dim][i] < pivot) {
            left_idx[lc]    = indexes[i];
            left_pts[0][lc] = points[0][i];
            left_pts[1][lc] = points[1][i];
            left_pts[2][lc] = points[2][i];
            lc += 1;
        } else {
            rght_idx[rc]    = indexes[i];
            rght_pts[0][rc] = points[0][i];
            rght_pts[1][rc] = points[1][i];
            rght_pts[2][rc] = points[2][i];
            rc += 1;
        }
    }

    length[0] = lc;
    length[1] = rc;
}

double classifyWithSlide(double **points, int *indexes,
                  int *left_idx, int *right_idx,
                  double **left_pts, double **rght_pts,
                  double pivot, int dim, int *length, int n)
/* Classify the points and slide the pivot if necessary*/
{
    classify(points, indexes, left_idx, right_idx, left_pts, rght_pts,
             pivot, dim, length, n);

    // Slide because empty left!
    if (length[0] == 0) {
        pivot = mind(points[dim], n) - 1e-10;
        classify(points, indexes, left_idx, right_idx, left_pts, rght_pts,
                 pivot, dim, length, n);

    // Slide because empty right!
    } else if (length[1] == 0) {
        pivot = maxd(points[dim], n) + 1e-10;
        classify(points, indexes, left_idx, right_idx, left_pts, rght_pts,
                 pivot, dim, length, n);
    }

    return pivot;
}
