#include <stdio.h>
#include <stdlib.h>
#include "io.h"
#include "definitions.h"

// Read libraries and dependencies for UNIX
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
//

int STmxyz_coarse()
{
  // Allocate the arrays in memory
  double (*mxyz_coarse)[NSD] = (double(*)[NSD]) malloc(sizeof(double)*nsd*nn_coarse);
  read_nodes(&mxyz_coarse[0][0], nn_coarse, nsd, "Data/mxyz.coarse");

  int i, j;

  // mxyz_coarse
  double **mxyz_coarseT = malloc(nsd*sizeof(double*));
  for (i=0; i<nsd; i++) {
    mxyz_coarseT[i] = malloc(nn_coarse*sizeof(double));
  }
  for (i=0; i<nn_coarse; i++) {
    for (j=0; j<nsd; j++) {
      mxyz_coarseT[j][i] = mxyz_coarse[i][j];
    }
  }


  FILE *f;
  f = fopen("Data/mxyz_ns.coarse", "wb");
  for (int i=0; i<nsd; i++) {
    fwrite( mxyz_coarseT[i], nn_coarse, sizeof(double), f );
  }
  fclose(f);

  return 0;
}





int STmien_coarse()
{
  int    (*mien_coarse)[NEN] = (int(*)[NEN])    malloc(sizeof(int)*nen*ne_coarse);
  read_connectivity(&mien_coarse[0][0], ne_coarse, nen, "Data/mien.coarse");

  int i, j;

  // mien_coarse
  int **mien_coarseT = malloc(nen*sizeof(int*));
  for (i=0; i<nen; i++) {
    mien_coarseT[i] = malloc(ne_coarse*sizeof(int));
  }
  for (i=0; i<ne_coarse; i++) {
    for (j=0; j<nen; j++) {
      mien_coarseT[j][i] = mien_coarse[i][j];
    }
  }

  FILE *f;
  f = fopen("Data/mien_ns.coarse", "wb");
  for (int i=0; i<nen; i++) {
    fwrite( mien_coarseT[i], ne_coarse, sizeof(int), f );
  }
  fclose(f);

  return 0;
}





int STdata_coarse()
{
  double (*data_coarse)[NDF] = (double(*)[NDF]) malloc(sizeof(double)*ndf*nn_coarse);
  read_data(        &data_coarse[0][0], nn_coarse, ndf, "Data/pres.coarse");

  int i, j;

  // data_coarse
  double **data_coarseT = malloc(ndf*sizeof(double*));
  for (i=0; i<ndf; i++) {
    data_coarseT[i] = malloc(nn_coarse*sizeof(double));
  }
  for (i=0; i<nn_coarse; i++) {
    for (j=0; j<ndf; j++) {
      data_coarseT[j][i] = data_coarse[i][j];
    }
  }

  FILE *f;
  f = fopen("Data/data_ns.coarse", "wb");
  for (i=0; i<ndf; i++) {
    fwrite( data_coarseT[i], nn_coarse, sizeof(double), f );
  }
  fclose(f);

  return 0;
}


int STmxyz_fine()
{
  double (*mxyz_fine)[NSD] = (double(*)[nsd]) malloc(sizeof(double)*nsd*nn_fine);
  read_nodes(&mxyz_fine[0][0],   nn_fine,   nsd, "Data/mxyz.fine");

  int i, j;

  // mxyz_fine
  double **mxyz_fineT = malloc(nsd*sizeof(double*));
  for (i=0; i<nsd; i++) {
    mxyz_fineT[i] = malloc(nn_fine*sizeof(double));
  }
  for (i=0; i<nn_fine; i++) {
    for (j=0; j<nsd; j++) {
      mxyz_fineT[j][i] = mxyz_fine[i][j];
    }
  }


  FILE *f;
  f = fopen("Data/mxyz_ns.fine", "wb");
  for (int i=0; i<nsd; i++) {
    fwrite( mxyz_fineT[i], nn_fine, sizeof(double), f );
  }
  fclose(f);

  return 0;
}



int main()
{

  // printf("1\n");
  STmxyz_coarse();

  // printf("2\n");
  STmien_coarse();

  // printf("3\n");
  STdata_coarse();

  // printf("4\n");
  STmxyz_fine();

  return 0;
}
