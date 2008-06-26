/*  Circular queue stencil code
 *  Kaushik Datta (kdatta@cs.berkeley.edu)
 *  University of California Berkeley
 *
 *  This code implements the circular queue algorithm for stencil codes.
 *  Intermediate queues, each with three revolving planes, store temporary
 *  results until the final result is written to the target array.  Unlike
 *  the time skewing algorithm, this algorithm will perform redundant
 *  computation between adjacent slabs.
 *
 *  NOTE: Only the cache block's y-dimension is used in this code; it
 *  specifies the size of the circular queue's y-dimension.  The grid's
 *  y-dimension needs to be a multiple of the cache block's y-dimension.
 */

#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#define MAX(x,y) (x > y ? x : y)

double *queuePlanes, *queuePlane0, *queuePlane1, *queuePlane2;
int *queuePlanesIndices;

/* This method creates the circular queues that will be needed for the
   circular_queue() method.  It is only called when more than one iteration
   is being performed. */
void CircularQueueInit(int nx, int ty, int timesteps) {
  int numPointsInQueuePlane, t;
  
  queuePlanesIndices = (int *) malloc((timesteps-1) * sizeof(int));
  
  if (queuePlanesIndices==NULL) {
    printf("Error on array queuePlanesIndices malloc.\n");
    exit(EXIT_FAILURE);
  }
  
  int queuePlanesIndexPtr = 0;
  
  for (t=1; t < timesteps; t++) {
    queuePlanesIndices[t-1] = queuePlanesIndexPtr;
    numPointsInQueuePlane = (ty+2*(timesteps-t)) * nx;
    queuePlanesIndexPtr += numPointsInQueuePlane;
  }

  queuePlanes = (double *) malloc(3 * queuePlanesIndexPtr * sizeof(double));
  
  if (queuePlanes==NULL) {
    printf("Error on array queuePlanes malloc.\n");
    exit(EXIT_FAILURE);
  }
  
  queuePlane0 = queuePlanes;
  queuePlane1 = &queuePlanes[queuePlanesIndexPtr];
  queuePlane2 = &queuePlanes[2 * queuePlanesIndexPtr];
}

/* This method traverses each slab and uses the circular queues to perform the
   specified number of iterations.  The circular queue at a given timestep is
   shrunken in the y-dimension from the circular queue at the previous timestep. */
#ifdef STENCILTEST
void StencilProbe_circqueue(double *A0, double *Anext, int nx, int ny, int nz,
  int tx, int ty, int tz, int timesteps) {
#else
void StencilProbe(double *A0, double *Anext, int nx, int ny, int nz,
  int tx, int ty, int tz, int timesteps) {
#endif
  double *readQueuePlane0, *readQueuePlane1, *readQueuePlane2, *writeQueuePlane, *tempQueuePlane;
  int blockMin_y, blockMax_y;
  int writeBlockMin_y, writeBlockMax_y;
  int writeBlockRealMin_y, writeBlockRealMax_y;
  int readBlockUnitStride_y, writeBlockUnitStride_y;
  int readOffset, writeOffset;
  int i, j, k, s, t;
  
  double fac = A0[0];
  int numBlocks_y = (ny-2)/ty;

  for (s=0; s < numBlocks_y; s++) {
    for (k=1; k < (nz+timesteps-2); k++) {
      for (t=0; t < timesteps; t++) {
	if ((k > t) && (k < (nz+t-1))) {

	  if (t == 0) {
	    readQueuePlane0 = &A0[Index3D(nx, ny, 0, 0, k-1)];
	    readQueuePlane1 = &A0[Index3D(nx, ny, 0, 0, k)];
	    readQueuePlane2 = &A0[Index3D(nx, ny, 0, 0, k+1)];
	  }
	  else {
	    readQueuePlane0 = &queuePlane0[queuePlanesIndices[t-1]];
	    readQueuePlane1 = &queuePlane1[queuePlanesIndices[t-1]];
	    readQueuePlane2 = &queuePlane2[queuePlanesIndices[t-1]];
	  }

	  // determine the edges of the queues
	  writeBlockMin_y = s * ty - (timesteps-t) + 2;
	  writeBlockMax_y = (s+1) * ty + (timesteps-t);
	  writeBlockRealMin_y = writeBlockMin_y;
	  writeBlockRealMax_y = writeBlockMax_y;

	  if (writeBlockMin_y < 1) {
	    writeBlockMin_y = 0;
	    writeBlockRealMin_y = 1;
	  }
	  if (writeBlockMax_y > (ny-1)) {
	    writeBlockMax_y = ny;
	    writeBlockRealMax_y = ny-1;
	  }

	  if (t == (timesteps-1)) {
	    writeQueuePlane = Anext;
	    writeOffset = 0;
	  }
	  else {
	    writeQueuePlane = &queuePlane2[queuePlanesIndices[t]];
	    writeOffset = Index3D(nx, ny, 0, writeBlockMin_y, k-t);
	  }

	  if ((writeBlockMin_y == 0) || (t == 0)) {
	    readOffset = Index3D(nx, ny, 0, 0, k-t);
	  }
	  else {
	    readOffset = Index3D(nx, ny, 0, writeBlockMin_y-1, k-t);
	  }

	  // use ghost cells for the bottommost and topmost planes
	  if (k == (t+1)) {
	    readQueuePlane0 = A0;
	  }
	  if (k == (nz+t-2)) {
	    readQueuePlane2 = &A0[Index3D(nx, ny, 0, 0, nz-1)];
	  }

	  // copy ghost cells
	  if (t < (timesteps-1)) {
	    for (j=(writeBlockMin_y+1); j < (writeBlockMax_y-1); j++) {
	      writeQueuePlane[Index3D(nx, ny, 0, j, k-t) - writeOffset] = readQueuePlane1[Index3D(nx, ny, 0, j, k-t) - readOffset];
	      writeQueuePlane[Index3D(nx, ny, nx-1, j, k-t) - writeOffset] = readQueuePlane1[Index3D(nx, ny, nx-1, j, k-t) - readOffset];
	    }
	    if (writeBlockMin_y == 0) {
	      for (i=1; i < (nx-1); i++) {
		writeQueuePlane[Index3D(nx, ny, i, writeBlockMin_y, k-t) - writeOffset] = readQueuePlane1[Index3D(nx, ny, i, writeBlockMin_y, k-t) - readOffset];
	      }
	    }
	    if (writeBlockMax_y == ny) {
	      for (i=1; i < (nx-1); i++) {
		writeQueuePlane[Index3D(nx, ny, i, writeBlockRealMax_y, k-t) - writeOffset] = readQueuePlane1[Index3D(nx, ny, i, writeBlockRealMax_y, k-t) - readOffset];
	      }
	    }
	  }

	  // actual calculations
	  for (j=writeBlockRealMin_y; j < writeBlockRealMax_y; j++) {
	    for (i=1; i < (nx-1); i++) {
	      writeQueuePlane[Index3D(nx, ny, i, j, k-t) - writeOffset] = 
		readQueuePlane0[Index3D(nx, ny, i, j, k-t) - readOffset] +
		readQueuePlane2[Index3D(nx, ny, i, j, k-t) - readOffset] +
		readQueuePlane1[Index3D(nx, ny, i, j-1, k-t) - readOffset] +
		readQueuePlane1[Index3D(nx, ny, i-1, j, k-t) - readOffset] +
		readQueuePlane1[Index3D(nx, ny, i+1, j, k-t) - readOffset] +
		readQueuePlane1[Index3D(nx, ny, i, j+1, k-t) - readOffset]
		- 6.0 * readQueuePlane1[Index3D(nx, ny, i, j, k-t) - readOffset] / (fac*fac);
	    }
	  }
	}
      }
      if (t > 0) {
	tempQueuePlane = queuePlane0;
	queuePlane0 = queuePlane1;
	queuePlane1 = queuePlane2;
	queuePlane2 = tempQueuePlane;
      }
    }
  }
}
