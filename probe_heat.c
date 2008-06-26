/*
	StencilProbe Heat Equation
	Implements 7pt stencil from Chombo's heattut example.
*/
#include <stdio.h>
#include "common.h"

#ifdef STENCILTEST
void StencilProbe_naive(double* A0, double* Anext, int nx, int ny, int nz,
			int tx, int ty, int tz, int timesteps) {
#else
void StencilProbe(double *A0, double *Anext, int nx, int ny, int nz,
                  int tx, int ty, int tz, int timesteps) {
#endif
  // Fool compiler so it doesn't insert a constant here
  double fac = A0[0];
  double *temp_ptr;
  int i, j, k, t;
  
  for (t = 0; t < timesteps; t++) {
    for (k = 1; k < nz - 1; k++) {
      for (j = 1; j < ny - 1; j++) {
	for (i = 1; i < nx - 1; i++) {
	  Anext[Index3D (nx, ny, i, j, k)] = 
	    A0[Index3D (nx, ny, i, j, k + 1)] +
	    A0[Index3D (nx, ny, i, j, k - 1)] +
	    A0[Index3D (nx, ny, i, j + 1, k)] +
	    A0[Index3D (nx, ny, i, j - 1, k)] +
	    A0[Index3D (nx, ny, i + 1, j, k)] +
	    A0[Index3D (nx, ny, i - 1, j, k)]
	    - 6.0 * A0[Index3D (nx, ny, i, j, k)] / (fac*fac);
	}
      }
    }
    temp_ptr = A0;
    A0 = Anext;
    Anext = temp_ptr;
  }
}
