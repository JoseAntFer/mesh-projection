#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
// #include "interfaces.h"

// Read libraries and dependencies for UNIX
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
//

// Use this routine in order to swap bytes
void swapbytes(char* array, int nelem, int elemsize) {
  register int sizet, sizem, i, j;
  char* bytea;
  sizet = elemsize;
  sizem = sizet - 1;
  bytea = malloc(sizet);
  for(i = 0; i < nelem; ++i) {
    memcpy((void*) bytea, (void*) (array + i*sizet), sizet);
    for(j = 0; j < sizet; ++j) {
      array[i*sizet + j] = bytea[sizem - j];
    }
  }
  free(bytea);
}

// Use this routine in order to read a file containing coordinates (mxyz)
void read_nodes(
    double* mxyz,     // Output: An array containing the coordinates
    unsigned nn,      // Input: The number of nodes to read
    unsigned nsd,     // Input: The number of space dimensions
    char* file_name   // Input: The path to the coordinate file
    ){

    int  f = open(file_name, O_RDONLY);
    if (!f) { printf("Unable to open file!"); }

    read(f, mxyz, 8*nn*nsd);
    swapbytes( (char *)mxyz, nn*nsd, sizeof(double) );
    close(f);
}


#define NEN 4
#define NDF 1
#define NSD 3
const unsigned nen = 4;
const unsigned ndf = 1;
const unsigned nsd = 3;
const unsigned nn_fine = 110618;
const unsigned nn_coarse = 26148;
const unsigned ne_coarse = 129968;



double **transpose(double **A, int m, int n)
{
  int i, j;
  double **B = malloc(n*sizeof(double*));

  for (i=0; i<n; i++) {B[i] = malloc(m*sizeof(double));}

  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      B[j][i] = A[i][j];
    }
  }

  return B;
}


int main()
{

  int i, j;
  double (*mxyz_fine)[NSD] = malloc(sizeof(double)*NSD*nn_fine);
  read_nodes(&mxyz_fine[0][0], nn_fine, NSD, "Data/mxyz.fine");

  for (i=0; i<10; i++) {
    for (j=0; j<NSD; j++) {
      printf("%.3f, ", mxyz_fine[i][j]);
    }
    printf("\n");
  }

  for (i=0; i<10; i++) { printf("%.3f, ", *mxyz_fine[i]); }


  // double **B = transpose(points, nn_fine, NSD);
  //
  // for (i=0; i<NSD; i++) {
  //   for (j=0; j<10; j++) {
  //     printf("%.3f, ", B[i][j]);
  //   }
  //   printf("\n");
  // }


  return 0;
}
