/*
	Stencil Probe
	Main function.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "util.h"
#include "cycle.h"
#ifdef HAVE_PAPI
#include <papi.h>
#endif
/* run.h has the run parameters */
#include "run.h"

void StencilProbe(double* A0, double* Anext, int nx, int ny, int nz,
                  int tx, int ty, int tz, int timesteps);

int main(int argc,char *argv[])
{
   double *Anext_naive, *Anext_test;
   double *A0_naive, *A0_test;
   int nx,ny,nz,tx,ty,tz,timesteps;
   int i,j;

	ticks t1, t2;
	double spt;

  /* parse command line options */
  if (argc < 8)
  {
    printf("Usage: %s <grid x> <grid y> <grid z> <block x> <block y> <block z> <timesteps>\n", argv[0]);
    return 1;
  }
  
  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  nz = atoi(argv[3]);
  tx = atoi(argv[4]);
  ty = atoi(argv[5]);
  tz = atoi(argv[6]);
  timesteps = atoi(argv[7]);
  printf("%dx%dx%d, blocking: %dx%dx%d, timesteps: %d\n",
    nx,ny,nz,tx,ty,tz,timesteps);

#ifdef HAVE_PAPI
	PAPI_library_init(PAPI_VER_CURRENT);
#endif

	/* find conversion factor from ticks to seconds */
	spt = seconds_per_tick();

  	/* allocate arrays */ 
   Anext_naive=(double*)malloc(sizeof(double)*nx*ny*nz);
   A0_naive=(double*)malloc(sizeof(double)*nx*ny*nz);
   Anext_test=(double*)malloc(sizeof(double)*nx*ny*nz);
   A0_test=(double*)malloc(sizeof(double)*nx*ny*nz);

	printf("USING TIMER: %s \t  SECONDS PER TICK:%g \n", TIMER_DESC, spt);
	
	/* initialize arrays to all ones */
   	StencilInit(nx,ny,nz,Anext_naive);
   	StencilInit(nx,ny,nz,A0_naive);
  /* copy into test arrays */
    for (j=0; j<nx*ny*nz; j++)
    {
      Anext_test[j] = Anext_naive[j];
      A0_test[j] = A0_naive[j];
    }
	/* stencil function */ 
  printf("Checking Rivera (single-timestep) blocking...\n");
	StencilProbe_naive(A0_naive, Anext_naive, nx, ny, nz, tx, ty, tz, timesteps);
  StencilProbe_rivera(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps);
  check_vals(Anext_naive, Anext_test, nx, ny, nz);

	/* initialize arrays to all ones */
   	StencilInit(nx,ny,nz,Anext_naive);
   	StencilInit(nx,ny,nz,A0_naive);
  /* copy into test arrays */
    for (j=0; j<nx*ny*nz; j++)
    {
      Anext_test[j] = Anext_naive[j];
      A0_test[j] = A0_naive[j];
    }
	/* stencil function */ 
  printf("Checking Cache-Oblivious blocking...\n");
	StencilProbe_naive(A0_naive, Anext_naive, nx, ny, nz, tx, ty, tz, timesteps);
  StencilProbe_oblivious(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps);
  check_vals(Anext_naive, Anext_test, nx, ny, nz);



	/* free arrays */
   free(Anext_naive);
   free(A0_naive);
}

void check_vals(double* A, double* B, int nx, int ny, int nz)
{
  int i, same, different;
  same=different=0;

  for (i=0; i<nx*ny*nz; i++)
    if (A[i] == B[i] || fabs(A[i]-B[i]) < 0.001  )
      same++;
    else
      different++;

  printf("Same: %d   Different: %d\n", same, different);

}
