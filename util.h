#ifndef _PROBE_H_
#define _PROBE_H_

#include "common.h"


/*
  This initializes the array A to be all 1's.  
  This is nearly superfluous (could use memset), but
  provides convenience and consistency nonetheless...
 */
void StencilInit(int nx,int ny,int nz, /* size of the array */
		 double *A); /* the array to initialize to 1's */


void clear_cache();

double seconds_per_tick();

#endif
