/*
	Stencil Probe
	Main function.
*/

#include <stdio.h>
#include <stdlib.h>
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
  double *Anext;
  double *A0;
  int nx,ny,nz,tx,ty,tz,timesteps;
  int i;
  
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
  
  /* allocate arrays */ 
  Anext=(double*)malloc(sizeof(double)*nx*ny*nz);
  A0=(double*)malloc(sizeof(double)*nx*ny*nz);
  
  printf("USING TIMER: %s \t  SECONDS PER TICK:%g \n", TIMER_DESC, spt);
  
  for (i=0;i<NUM_TRIALS;i++) {
    /* initialize arrays to all ones */
    StencilInit(nx,ny,nz,Anext);
    StencilInit(nx,ny,nz,A0);
#ifdef CIRCULARQUEUEPROBE
    if (timesteps > 1) {                                                                                                                      
      CircularQueueInit(nx, ty, timesteps);                                                                                                   
    }
#endif    

    // clear_cache();
    
    t1 = getticks();	
    
    /* stencil function */ 
    StencilProbe(A0, Anext, nx, ny, nz, tx, ty, tz, timesteps);
    
    t2 = getticks();
    
    printf("elapsed ticks: %g  time:%g \n", elapsed(t2, t1), spt * elapsed(t2,t1));
  }
  
  /* free arrays */
  free(Anext);
  free(A0);
  return EXIT_SUCCESS;
}
