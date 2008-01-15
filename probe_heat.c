/*
	StencilProbe Heat Equation
	Impelements 7pt stencil from Chombo's heattut example.
*/
#include "common.h"

void StencilProbe(double* A0, double* Anext, int nx, int ny, int nz,
                  int tx, int ty, int tz, int timesteps)
{
  /* Fool compiler so it doesn't insert a constant here */
  double fac = A0[0];
  double *myA0, *myAnext;
  int t, i, j, k;

  for (t = 0; t < timesteps; t++)
  {
    if (!(t%2))
    { myA0 = A0; myAnext = Anext; }
    else
    { myA0 = Anext; myAnext = A0; }
    /* Standard 3 nested loops stencil */
   for (k = 1; k < nz - 1; k++)
    {
      for (j = 1; j < ny - 1; j++)
	{
	  for (i = 1; i < nx - 1; i++)
	    {
		myAnext[Index3D (nx, ny, i, j, k)] = 
			myA0[Index3D (nx, ny, i, j, k + 1)] +
		    myA0[Index3D (nx, ny, i, j, k - 1)] +
		    myA0[Index3D (nx, ny, i, j + 1, k)] +
		    myA0[Index3D (nx, ny, i, j - 1, k)] +
		    myA0[Index3D (nx, ny, i + 1, j, k)] +
		    myA0[Index3D (nx, ny, i - 1, j, k)]
			- 6.0 * myA0[Index3D (nx, ny, i, j, k)] / (fac*fac);
	    }
	}
    }
    }
}
