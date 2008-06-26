/*  Time skewing stencil code
 *  Kaushik Datta (kdatta@cs.berkeley.edu)
 *  University of California Berkeley
 *
 *  This code implements the time skewing method.  The cache blocks need to be
 *  traversed in a specific order for the algorithm to work properly.
 *
 *  NOTE: The number of iterations can only be up to one greater than the
 *  smallest cache block dimension.  If you wish to do more iterations, there
 *  are two options:
 *    1.  Make the smallest cache block dimension larger.
 *    2.  Split the number of iterations into smaller runs where each run
 *        conforms to the above rule.
 */
#include "common.h"
#define MAX(x,y) (x > y ? x : y)

/* This method traverses all of the cache blocks in a specific order to preserve
   dependencies.  For each cache block, it performs (possibly) several iterations while
   still respecting boundary conditions.
   NOTE: Positive slopes indicate that each iteration goes further out from the center
   of the current cache block, while negative slopes go toward the block center. */
#ifdef STENCILTEST
void StencilProbe_timeskew(double *A0, double *Anext, int nx, int ny, int nz,
  int tx, int ty, int tz, int timesteps) {
#else
void StencilProbe(double *A0, double *Anext, int nx, int ny, int nz,
  int tx, int ty, int tz, int timesteps) {
#endif
  double fac = A0[0];
  double *temp_ptr;
  double *myA0, *myAnext;

  int neg_x_slope, pos_x_slope, neg_y_slope, pos_y_slope, neg_z_slope, pos_z_slope;
  int blockMin_x, blockMin_y, blockMin_z;
  int blockMax_x, blockMax_y, blockMax_z;
  int ii, jj, kk, i, j, k, t;

  for (kk=1; kk < nz-1; kk+=tz) {
    neg_z_slope = 1;
    pos_z_slope = -1;

    if (kk == 1) {
      neg_z_slope = 0;
    }
    if (kk == nz-tz-1) {
      pos_z_slope = 0;
    }
    for (jj=1; jj < ny-1; jj+=ty) {
      neg_y_slope = 1;
      pos_y_slope = -1;
      
      if (jj == 1) {
	neg_y_slope = 0;
      }
      if (jj == ny-ty-1) {
	pos_y_slope = 0;
      }
      for (ii=1; ii < nx-1; ii+=tx) {
	neg_x_slope = 1;
	pos_x_slope = -1;
	
	if (ii == 1) {
	  neg_x_slope = 0;
	}
	if (ii == nx-tx-1) {
	  pos_x_slope = 0;
	}

	myA0 = A0;
	myAnext = Anext;
	
	for (t=0; t < timesteps; t++) {
	  blockMin_x = MAX(1, ii - t * neg_x_slope);
	  blockMin_y = MAX(1, jj - t * neg_y_slope);
	  blockMin_z = MAX(1, kk - t * neg_z_slope);
	  
	  blockMax_x = MAX(1, ii + tx + t * pos_x_slope);
	  blockMax_y = MAX(1, jj + ty + t * pos_y_slope);
	  blockMax_z = MAX(1, kk + tz + t * pos_z_slope);
	  
	  for (k=blockMin_z; k < blockMax_z; k++) {
	    for (j=blockMin_y; j < blockMax_y; j++) {
	      for (i=blockMin_x; i < blockMax_x; i++) {
		myAnext[Index3D (nx, ny, i, j, k)] = 
		  myA0[Index3D (nx, ny, i, j, k+1)] +
		  myA0[Index3D (nx, ny, i, j, k-1)] +
		  myA0[Index3D (nx, ny, i, j+1, k)] +
		  myA0[Index3D (nx, ny, i, j-1, k)] +
		  myA0[Index3D (nx, ny, i+1, j, k)] +
		  myA0[Index3D (nx, ny, i-1, j, k)]
		  - 6.0 * myA0[Index3D (nx, ny, i, j, k)] / (fac*fac);
	      }
	    }
	  }
	  temp_ptr = myA0;
	  myA0 = myAnext;
	  myAnext = temp_ptr;
	}
      }
    }
  }
}
