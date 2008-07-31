/*
	Stencil Probe utilities
	Helper functions for the probe.
*/
#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "cycle.h"



/*
  This initializes the array A to be all 1's.  
  This is nearly superfluous (could use memset), but
  provides convenience and consistency nonetheless...
 */
void StencilInit(int nx,int ny,int nz, /* size of the array */
		 double *A){ /* the array to initialize to 1s */
  long last = nx*ny*nz;
  long i;

  for(i=0;i<last;i++) 
#ifdef RANDOMVALUES
	A[i]=(float)rand()/RAND_MAX;
#else
	A[i]=1.0;
#endif
}

/*
  This function determines ticks per second.
  Inspired by OSKI function (bebop.cs.berkeley.edu)
*/
double seconds_per_tick()
{
	ticks t0,t1;
	unsigned int i = 3;
	double spt = 0;

	while (spt <= 0)
	{

	t0=getticks();	
	sleep(i);
	t1=getticks();
	spt = (double)i / elapsed(t1,t0);
	i++;

	}

	return spt;
}

/*
  Function to clear the cache, preventing data items in cache
  from making subsequent trials run faster.
*/
void clear_cache()
{
  int i;
  float* tarray, accum;

  tarray = (float*) malloc(sizeof(float)*1310720);
  for (i=0,accum=0; i<1310719; i++)
    tarray[i] = 1.0;

}
