#ifndef KDTOOLS_H_INCLUDED
#define KDTOOLS_H_INCLUDED

double maxd(double *array, int n);
double mind(double *array, int n);
void boundingBox(double **points, double *bb1, double *bb2, int n, double tol);
int longestDimension(double *bb1, double *bb2);
void classify(double **points, int *indexes,
              int *left_idx, int *rght_idx,
              double **left_pts, double **rght_pts,
              double pivot, int dim, int *length, int n);
double classifyWithSlide(double **points, int *indexes,
              int *left_idx, int *right_idx,
              double **left_pts, double **rght_pts,
              double pivot, int dim, int *length, int n);



#endif // KDTOOLS_H_INCLUDED
