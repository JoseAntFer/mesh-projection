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

    int r = read(f, mxyz, 8*nn*nsd);
    swapbytes( (char *)mxyz, nn*nsd, sizeof(double) );
    close(f);
}

// Use this routine in order to read a file containing data (pres)
void read_data(
    double* data,     // Output: An array containing the data
    unsigned nn,      // Input: The number of nodes to read
    unsigned ndf,     // Input: The number of degrees of freedom
    char* file_name   // Input: The path to the data file
    ){

      int  f = open(file_name, O_RDONLY);
      if (!f) { printf("Unable to open file!"); }

      int r = read(f, data, 8*nn*ndf);
      swapbytes( (char *)data, nn*ndf, sizeof(double) );
      close(f);
}

// Use this routine in order to read the connectivity file (mien)
void read_connectivity(
    int* mien,        // Output: An array containing the connectivity information
    unsigned ne,      // Input: The number of elements
    unsigned nen,     // Input: The number of nodes per element
    char* file_name   // Input: The path to the connectivity file
    ){

    	int  f = open(file_name, O_RDONLY);
    	if (!f) { printf("Unable to open file!"); }

      int r = read(f, mien, 4*ne*nen);
      swapbytes( (char *)mien, ne*nen, sizeof(int) );
    	close(f);
}

// Use this routine in order to write the projected data (pres)
void write_data(
    double* data,     // Input: Data to be written
    unsigned nn,      // Input: The number of nodes to write
    unsigned ndf,     // Input: The number of degrees of freedom
    char* file_name   // Input: The path to the data file
    ){

      FILE *f;
    	f = fopen(file_name, "wb");
    	if (!f) { printf("Unable to open file!"); }

      swapbytes( (char *)data, nn*ndf, sizeof(double) );
    	fwrite( data, nn*ndf, sizeof(double), f );
    	fclose(f);
}

// Use this routine in order to read a file containing coordinates (mxyz)
double **read_d2D(
    unsigned d1,      // Input: The number of nodes to read
    unsigned d2,      // Input: The number of space dimensions
    char* file_name   // Input: The path to the coordinate file
    ){

    // Allocate ouput array
    int i, r;
    double **A = malloc(d1*sizeof(double*));
    for (i=0; i<d1; i++) {A[i] = malloc(d2*sizeof(double));}

    // Open file with data
    int f = open(file_name, O_RDONLY);
    if (!f) { printf("Unable to open file!"); }

    // Write data into the array
    for (i=0; i<d1; i++) {
      r = read(f, A[i], 8*d2);
    }
    close(f);

    return A;
}


// Use this routine in order to read a file containing coordinates (mxyz)
int **read_i2D(
    unsigned d1,      // Input: The number of nodes to read
    unsigned d2,      // Input: The number of space dimensions
    char* file_name   // Input: The path to the coordinate file
    ){

    // Allocate ouput array
    int i, r;
    int **A = malloc(d1*sizeof(int*));
    for (i=0; i<d1; i++) {A[i] = malloc(d2*sizeof(int));}

    // Open file with data
    int f = open(file_name, O_RDONLY);
    if (!f) { printf("Unable to open file!"); }

    // Write data into the array
    for (i=0; i<d1; i++) {
      r = read(f, A[i], 4*d2);
    }
    close(f);

    return A;
}
