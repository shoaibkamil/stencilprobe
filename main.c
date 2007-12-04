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
   const int nx=SIZE,ny=SIZE,nz=SIZE;
   int i;

	ticks t1, t2;
	double spt;

#ifdef HAVE_PAPI
	PAPI_library_init(PAPI_VER_CURRENT);
#endif

	/* find conversion factor from ticks to seconds */
	spt = seconds_per_tick();

  	/* allocate arrays */ 
   Anext=(double*)malloc(sizeof(double)*nx*ny*nz);
   A0=(double*)malloc(sizeof(double)*nx*ny*nz);

	printf("USING TIMER: %s \t  SECONDS PER TICK:%g \n", TIMER_DESC, spt);
	
   for (i=0;i<NUM_TRIALS;i++)
   {
	/* initialize arrays to all ones */
   	StencilInit(nx,ny,nz,Anext);
   	StencilInit(nx,ny,nz,A0);


	t1 = getticks();	
  	
	/* stencl function */ 
	StencilProbe(A0, Anext, nx, ny, nz, -1, -1, -1, 1);

	t2 = getticks();

	printf("elapsed ticks: %g  time:%g \n", elapsed(t2, t1), spt * elapsed(t2,t1));
   }

	/* free arrays */
   free(Anext);
   free(A0);
}
