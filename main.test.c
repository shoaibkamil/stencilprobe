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
void check_vals(double* A, double* B, int nx, int ny, int nz);

int main(int argc,char *argv[]) {
  double *A0_naive, *A0_test;
  double *Anext_naive, *Anext_test;
  double *Afinal_naive, *Afinal_test;
  int nx,ny,nz,tx,ty,tz,timesteps;
  int i,j;
  
  ticks t1, t2;
  double spt;
  
  /* parse command line options */
  if (argc < 8) {
    printf("\nUSAGE:\n%s <grid x> <grid y> <grid z> <block x> <block y> <block z> <timesteps>\n", argv[0]);
    printf("\nTIME SKEWING CONSTRAINTS:\nIn each dimension, <grid size - 2> should be a multiple of <block size>.\n");
    printf("\nCIRCULAR QUEUE CONSTRAINTS:\n<grid y - 2> should be a multiple of <block y>.  The block sizes in the other dimensions are ignored.\n\n");
    return EXIT_FAILURE;
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
  
  // allocate arrays
  A0_naive=(double*)malloc(sizeof(double)*nx*ny*nz);
  A0_test=(double*)malloc(sizeof(double)*nx*ny*nz);
  Anext_naive=(double*)malloc(sizeof(double)*nx*ny*nz);
  Anext_test=(double*)malloc(sizeof(double)*nx*ny*nz);

  // Run Naive Code
  StencilInit(nx,ny,nz,A0_naive);
  StencilInit(nx,ny,nz,Anext_naive);
  StencilProbe_naive(A0_naive, Anext_naive, nx, ny, nz, tx, ty, tz, timesteps);
  if (timesteps%2 == 0) {
    Afinal_naive = A0_naive;
  }
  else {
    Afinal_naive = Anext_naive;
  }
  
  printf("USING TIMER: %s \t  SECONDS PER TICK:%g \n", TIMER_DESC, spt);

  // Test Rivera Blocking
  StencilInit(nx,ny,nz,A0_test);
  StencilInit(nx,ny,nz,Anext_test);
  printf("Checking Rivera (single-timestep) blocking...\n");
  StencilProbe_rivera(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps);
  if (timesteps%2 == 0) {
    Afinal_test = A0_test;
  }
  else {
    Afinal_test = Anext_test;
  }
  check_vals(Afinal_naive, Afinal_test, nx, ny, nz);
  
  // Test Cache-Oblivious Blocking
  StencilInit(nx,ny,nz,A0_test);
  StencilInit(nx,ny,nz,Anext_test);
  printf("Checking Cache-Oblivious blocking...\n");
  StencilProbe_oblivious(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps);
  if (timesteps%2 == 0) {
    Afinal_test = A0_test;
  }
  else {
    Afinal_test = Anext_test;
  }
  check_vals(Afinal_naive, Afinal_test, nx, ny, nz);
  
  // Test Time-Skewed Blocking
  StencilInit(nx,ny,nz,A0_test);
  StencilInit(nx,ny,nz,Anext_test);
  printf("Checking Time-Skewed blocking...\n");
  StencilProbe_timeskew(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps);
  if (timesteps%2 == 0) {
    Afinal_test = A0_test;
  }
  else {
    Afinal_test = Anext_test;
  }
  check_vals(Afinal_naive, Afinal_test, nx, ny, nz);

  // Test Circular-Queue Blocking
  StencilInit(nx,ny,nz,A0_test);
  StencilInit(nx,ny,nz,Anext_test);
  if (timesteps > 1) {
    CircularQueueInit(nx, ty, timesteps);
  }
  printf("Checking Circular-Queue blocking...\n");
  StencilProbe_circqueue(A0_test, Anext_test, nx, ny, nz, tx, ty, tz, timesteps);
  Afinal_test = Anext_test;
  check_vals(Afinal_naive, Afinal_test, nx, ny, nz);
  
  /* free arrays */
  free(Anext_naive);
  free(A0_naive);
  free(Anext_test);
  free(A0_test);
  return EXIT_SUCCESS;
}

void check_vals(double* A, double* B, int nx, int ny, int nz) {
  int same, different;
  int i, j, k;

  same=different=0;

  for (k=0; k<nz; k++) {
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
	if (fabs(A[Index3D(nx,ny,i,j,k)] - B[Index3D(nx,ny,i,j,k)]) < 0.001)
	  same++;
	else {
	  different++;
	  printf("at index %d %d %d --- A: %3.2g, B %3.2g\n", i, j, k,
		 A[Index3D(nx,ny,i,j,k)], B[Index3D(nx,ny,i,j,k)]);
	}
      }
    }
  }
  printf("Same: %d   Different: %d\n", same, different);
}
