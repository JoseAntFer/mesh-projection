#ifndef IO_H_INCLUDED
#define IO_H_INCLUDED

void swapbytes(char* array, int nelem, int elemsize);

void read_nodes(
    double* mxyz,     // Output: An array containing the coordinates
    unsigned nn,      // Input: The number of nodes to read
    unsigned nsd,     // Input: The number of space dimensions
    char* file_name   // Input: The path to the coordinate file
);

void read_data(
   double* data,     // Output: An array containing the data
   unsigned nn,      // Input: The number of nodes to read
   unsigned ndf,     // Input: The number of degrees of freedom
   char* file_name   // Input: The path to the data file
);

void read_connectivity(
   int* mien,        // Output: An array containing the connectivity information
   unsigned ne,      // Input: The number of elements
   unsigned nen,     // Input: The number of nodes per element
   char* file_name   // Input: The path to the connectivity
);

void write_data(
   double* data,     // Input: Data to be written
   unsigned nn,      // Input: The number of nodes to write
   unsigned ndf,     // Input: The number of degrees of freedom
   char* file_name   // Input: The path to the data file
);

double **read_d2D(
  unsigned d1,      // Input: The number of nodes to read
  unsigned d2,     // Input: The number of space dimensions
  char* file_name   // Input: The path to the coordinate file
);

int **read_i2D(
    unsigned d1,      // Input: The number of nodes to read
    unsigned d2,      // Input: The number of space dimensions
    char* file_name   // Input: The path to the coordinate file
);

#endif // KDTREE_H_INCLUDED
