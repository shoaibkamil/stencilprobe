#define ds 1
#include "run.h"

#define idx(_i,_j,_k,_nx,_ny,_nz) ((_i)+(_nx)*((_j)+(_ny)*(_k)))

void walk3(double* A[], int nx, int ny, int nz,
           int t0, int t1, int x0, int dx0, int x1, int dx1,
           int y0, int dy0, int y1, int dy1,
           int z0, int dz0, int z1, int dz1)
{
    int dt = t1-t0;

    if (dt == 1 || (x1-x0)*(y1-y0)*(z1-z0) < CUTOFF)
    {
        int x,y,z,t;
        double fac = A[0][0];
				for (t=t0;t<t1;t++)
					for (z=z0+(t-t0)*dz0;z<z1+(t-t0)*dz1;z++)
              for (y=y0+(t-t0)*dy0;y<y1+(t-t0)*dy1;y++)
                  for (x=x0+(t-t0)*dx0;x<x1+(t-t0)*dx1;x++)
                {
                    A[(t+1)%2][idx(x,y,z,nx+2,ny+2,nz+2)] =
                        (A[t%2][idx((x+1),y,z,nx+2,ny+2,nz+2)]
                         + A[t%2][idx((x-1),y,z,nx+2,ny+2,nz+2)]
                         + A[t%2][idx(x,(y+1),z,nx+2,ny+2,nz+2)]
                         + A[t%2][idx(x,(y-1),z,nx+2,ny+2,nz+2)]
                         + A[t%2][idx(x,y,(z+1),nx+2,ny+2,nz+2)]
                         + A[t%2][idx(x,y,(z-1),nx+2,ny+2,nz+2)]
                         - 6.0*A[t%2][idx(x,y,z,nx+2,ny+2,nz+2)]) / (fac*fac);
                }
    }
    else if (dt > 1)
    {
        if (2* (z1-z0) + (dz1-dz0) * dt >= 4 * ds * dt)
        {
            int zm = (2* (z0+z1) + (2*ds+dz0+dz1) * dt) / 4;
            walk3(A,nx,ny,nz,t0,t1,x0,dx0,x1,dx1,y0,dy0,y1,dy1,z0,dz0,zm,-ds);
            walk3(A,nx,ny,nz,t0,t1,x0,dx0,x1,dx1,y0,dy0,y1,dy1,zm,-ds,z1,dz1);
        }
        else if (2* (y1-y0) + (dy1-dy0) * dt >= 4 * ds * dt)
        {
            int ym = (2* (y0+y1) + (2*ds+dy0+dy1) * dt) / 4;
            walk3(A,nx,ny,nz,t0,t1,x0,dx0,x1,dx1,y0,dy0,ym,-ds,z0,dz0,z1,dz1);
            walk3(A,nx,ny,nz,t0,t1,x0,dx0,x1,dx1,ym,-ds,y1,dy1,z0,dz0,z1,dz1);
        }
/*              else if (2* (x1-x0) + (dx1-dx0) * dt >= 4 * ds * dt)
                {
                int xm = (2* (x0+x1) + (2*ds+dx0+dx1) * dt) / 4;
                walk3(A,nx,ny,nz,t0,t1,x0,dx0,xm,-ds,y0,dy0,y1,dy1,z0,dz0,z1,dz1);
                walk3(A,nx,ny,nz,t0,t1,xm,-ds,x1,dx1,y0,dy0,y1,dy1,z0,dz0,z1,dz1);
                }
*/
        else
        {
            int s = dt/2;
            walk3(A,nx,ny,nz,t0,t0+s,x0,dx0,x1,dx1,y0,dy0,y1,dy1,z0,dz0,z1,dz1);
            walk3(A,nx,ny,nz,t0+s,t1,x0+dx0*s,dx0,x1+dx1*s,dx1,y0+dy0*s,dy0,y1+dy1*s,dy1,
                  z0+dz0*s,dz0,z1+dz1*s,dz1);
        }
    }
}

void StencilProbe(double* A0, double* Anext, int nx, int ny, int nz,
                  int tx, int ty, int tz, int timesteps)
{
    double* A[2] = {A0, Anext};
    int i;

    walk3(A, nx-2, ny-2, nz-2,
          0, timesteps,
          1, 0, nx-1, 0,
          1, 0, ny-1, 0,
          1, 0, nz-1, 0);
}

void probe_naive(int* nx, int* ny, int* nz, double* A0, double* Anext)
{
    double* A[2] = {A0, Anext};
    int x,y,z;
    int numx,numy,numz;
    double fac = 1.0;
    int t=0;

    numx = *nx;
    numy = *ny;
    numz = *nz;

    for (z=1; z<numz-1; z++)
        for (y=1; y<numy-1; y++)
            for (x=1; x<numx-1; x++)
                A[(t+1)%2][idx(x,y,z,numx,numy,numz)] =
                    (A[t%2][idx((x+1),y,z,numx,numy,numz)]
                     + A[t%2][idx((x-1),y,z,numx,numy,numz)]
                     + A[t%2][idx(x,(y+1),z,numx,numy,numz)]
                     + A[t%2][idx(x,(y-1),z,numx,numy,numz)]
                     + A[t%2][idx(x,y,(z+1),numx,numy,numz)]
                     + A[t%2][idx(x,y,(z-1),numx,numy,numz)]
                     - 6.0*A[t%2][idx(x,y,z,numx,numy,numz)]) / (fac*fac);
} 

