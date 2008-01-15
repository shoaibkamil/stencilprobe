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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#define Index3D(_i,_j,_k) ((_i)+(tx)*((_j)+(ty)*(_k)))


/* This method traverses all of the cache blocks in a specific order to preserve
   dependencies.  For each cache block, it performs (possibly) several iterations while
   still respecting boundary conditions.
   NOTE: Positive slopes indicate that each iteration goes further out from the center
   of the current cache block, while negative slopes go toward the block center. */
void StencilProbe(double *A0, double *Anext, int nx, int ny, int nz,
  int tx, int ty, int tz, int timesteps) {

  /* Formerly globals */
  uint32_t numBlocks_x = nx / tx;
  uint32_t numBlocks_y = ny / ty;
  uint32_t numBlocks_z = nz / tz;
  uint32_t numBlockCells_x = tx;
  uint32_t numBlockCells_y = ty;
  uint32_t numBlockCells_z = tz;
  uint32_t realMin_x = 1;
  uint32_t realMin_y = 1;
  uint32_t realMin_z = 1;
  uint32_t realMax_x = realMin_x + nx;
  uint32_t realMax_y = realMin_y + ny;
  uint32_t realMax_z = realMin_z + nz;
  uint32_t numIter = timesteps;
  double *A = A0;
  double *B = Anext;
  double fac = A0[0];

  double * __restrict__ temp_ptr;
  uint32_t block_i, block_j, block_k;
  uint32_t neg_x_slope, pos_x_slope, neg_y_slope, pos_y_slope, neg_z_slope, pos_z_slope;
  uint32_t initBlockMin_x, initBlockMin_y, initBlockMin_z;
  uint32_t initBlockMax_x, initBlockMax_y, initBlockMax_z;
  uint32_t blockMin_x, blockMin_y, blockMin_z;
  uint32_t blockMax_x, blockMax_y, blockMax_z;
  uint32_t i, j, k, t;

  for (block_i=0; block_i < numBlocks_x; block_i++) {
    neg_x_slope = 1;
    pos_x_slope = -1;
    
    if (block_i == 0) {
      neg_x_slope = 0;
    }
    if (block_i == (numBlocks_x-1)) {
      pos_x_slope = 0;
    }
    
    for (block_j=0; block_j < numBlocks_y; block_j++) {
      neg_y_slope = 1;
      pos_y_slope = -1;
      
      if (block_j == 0) {
	neg_y_slope = 0;
      }
      if (block_j == (numBlocks_y-1)) {
	pos_y_slope = 0;
      }
      
      for (block_k=0; block_k < numBlocks_z; block_k++) {
	neg_z_slope = 1;
	pos_z_slope = -1;
	
	if (block_k == 0) {
	  neg_z_slope = 0;
	}
	if (block_k == (numBlocks_z-1)) {
	  pos_z_slope = 0;
	}
	
	double * __restrict__ Read_start = A;
	double * __restrict__ Write_start = B;
	
	initBlockMin_x = realMin_x + block_i * numBlockCells_x;
	initBlockMin_y = realMin_y + block_j * numBlockCells_y;
	initBlockMin_z = realMin_z + block_k * numBlockCells_z;
	
	initBlockMax_x = initBlockMin_x + numBlockCells_x;
	initBlockMax_y = initBlockMin_y + numBlockCells_y;
	initBlockMax_z = initBlockMin_z + numBlockCells_z;
	
	for (t=0; t < numIter; t++) {
	  // set block size based on number of iterations performed
	  blockMin_x = initBlockMin_x - t * neg_x_slope;
	  blockMin_y = initBlockMin_y - t * neg_y_slope;
	  blockMin_z = initBlockMin_z - t * neg_z_slope;
	  
	  blockMax_x = initBlockMax_x + t * pos_x_slope;
	  blockMax_y = initBlockMax_y + t * pos_y_slope;
	  blockMax_z = initBlockMax_z + t * pos_z_slope;

	  // actual calculations
	  for (i=blockMin_x; i < blockMax_x; i++) {
	    for (j=blockMin_y; j < blockMax_y; j++) {
	      for (k=blockMin_z; k < blockMax_z; k++) {
		Write_start[Index3D(i,j,k)] = -6.0 * Read_start[Index3D(i,j,k)] / fac*fac + 
		  (Read_start[Index3D(i-1,j,k)] + Read_start[Index3D(i+1,j,k)] 
      + Read_start[Index3D(i,j-1,k)] + Read_start[Index3D(i,j+1,k)] 
      + Read_start[Index3D(i,j,k-1)] + Read_start[Index3D(i,j,k+1)]);
	      }
	    }
	  }
	  temp_ptr = Read_start;
	  Read_start = Write_start;
	  Write_start = temp_ptr;
	}
      }
    }
  }
}

/*
int main(int argc, char *argv[]) {
  double results[NUM_TRIALS];
#if !defined(DEBUG)
#if defined(PAPI_ENABLED)
  int papi_setnum, num_desired, num_sets;
#else
  double median_counts_per_sec;
#endif
#endif

  printf("\n7-point stencil, no add, time skewed C code with non-periodic boundary conditions\n");

  // initialize arrays
  init_flush_cache_array();
  malloc_grids(argv);

#if !defined(DEBUG)
#if defined(PAPI_ENABLED)
  // initialize papi
  int desired_events[] = {PAPI_TOT_CYC, PAPI_FP_INS, PAPI_L2_DCA, PAPI_L2_DCM, PAPI_L3_DCM, PAPI_TLB_DM, PAPI_LD_INS, PAPI_SR_INS};
  num_desired = 9;
  PAPI_event_set_wrapper_t* event_sets;
  papi_init(desired_events, num_desired, &event_sets, &num_sets);
#else
  // calculate clock rate
  GET_CLOCK_RATE(results, NUM_TRIALS);
  median_counts_per_sec = find_median(results, NUM_TRIALS);
  printf("Median ticks per second = %e\n", median_counts_per_sec);
#endif
#endif

  printf("\n");

#if defined(DEBUG)
  init_grids();
  printf("Time skewing:\n");
  printf("\nGRID A BEFORE:");
  print_grid(A);
  printf("\nGRID B BEFORE:");
  print_grid(B);

  time_skewing();

  printf("\nGRID A AFTER:");
  print_grid(A);
  printf("\nGRID B AFTER:");
  print_grid(B);
#else
#if defined(PAPI_ENABLED)
  printf("Time skewing:\n");
  for (papi_setnum=0; papi_setnum < num_sets; papi_setnum++) {
    PAPI_MAKE_MEASUREMENTS(event_sets[papi_setnum].set, time_skewing(), NUM_TRIALS, results);
    print_papi_measurements(&(event_sets[papi_setnum]), results, NUM_TRIALS);
  }
  printf("\n");
#else
  printf("Time skewing:\n");
  TIMER_MAKE_MEASUREMENTS(time_skewing(), results, NUM_TRIALS);
  print_timer_measurements(results, NUM_TRIALS, median_counts_per_sec);
  printf("\n");
#endif
#endif

  printf("\nFinal interior values: A[%d, %d, %d] = %lf, B[%d, %d, %d] = %lf\n", nx/2, ny/2, nz/2, A[Index3D(nx/2, ny/2, nz/2)], nx/2, ny/2, nz/2, B[Index3D(nx/2, ny/2, nz/2)]);
  fc_checksum();
  free_grids();
  return EXIT_SUCCESS;
}*/
