#include <stdio.h>

#include <cublas.h>

#include "mp1-util.h"

#include <vector>


extern "C"{
#include "trajectory.h"
}

#define MAXBINS 500000
#define BINSIZE128 64
#define BINSIZE256 256

#ifndef PI
#define PI 3.141592653f
#endif

#define DBG_DCF 1
#define EPS 2.2204e-16

enum ConvType {ctCONST=0, ctSEP};

const char savePath[] = "/home_local/noaddy/research/data/code/dcfjcu/";
enum CompartmentSize{
  COMPARTMENT_64=64, 
  COMPARTMENT_128=128
};

Timer timer;

FILE* file;
char filename[128];

/**
 \param[in]	k	Distance from k-space origin
 \param[in]	F	k-space scaling factor
*/
__device__ float P2D(float k, float F)
{
	float val = 1.0f;
	
	if(k>1.0e-6)
	{
		val = F*j1f(F*PI*k)/(2.0f*k)/(PI*F*F/4.0f);
		val *= val;
	}
	
	return val;
}

/**
 \param[in]	k	Distance from k-space origin
 \param[in]	F	k-space scaling factor
*/
__device__ float P3D(float k, float F)
{
	float val;
// 	float s, c;
	float pfk;
	
	/*if(k<1e-6)
		val = 1.0f;
	else
	{*/
		pfk = PI*F*k;
// 		s = sinf(pfk);
// 		c = cosf(pfk);
// 		val = (s-pfk*c)/(2.0f*PI*PI*k*k*k) * (s-pfk*c)/(2.0f*PI*PI*k*k*k)/Pmax/Pmax;
		val = 3.0f*(sinf(pfk) - pfk*cosf(pfk))/(pfk*pfk*pfk);
		val *= val;
	//}
	
	return val;
}

__device__ float P3D1(float k, float F)
{
	
	return 9.0f*(sinf(PI*F*k) - PI*F*k*cosf(PI*F*k))*(sinf(PI*F*k) - PI*F*k*cosf(PI*F*k))/(PI*F*k*PI*F*k*PI*F*k * PI*F*k*PI*F*k*PI*F*k);
}

__device__ float P1D(float k, float F)
{
	float pfk = PI*F*k;
	float val;
	
// 	Avoid NaN errors
	if(k>EPS)
		val = sinf(pfk)/pfk;
	else
		val = 0.0f;

	return val*val;
}

/**
 \param[in]	connections	
 \param[in]	ndim	Number of dimensions	
*/
__global__ void conv_kernel(unsigned int *connections, int numConnections, int *bins, float *ks, int ndim, float kernelRad, float F, float* Win, float *Wout)
{
	unsigned int index = blockIdx.x + gridDim.x*blockIdx.y;
	unsigned int tx = threadIdx.x;
// 	unsigned int ty = threadIdx.y;
	unsigned int b1, b2;
	float3 k1, k2;
	float d, d2;
	float P;
	int n1, n2;
	int i2;
	
	if(index>=numConnections)
		return;
	
	b1 = connections[2*index];
	b2 = connections[2*index+1];
	
	n1 = bins[BINSIZE128*b1+tx];
	if(n1<0)
		return;
		
	if(b1==b2)
	{
		k1.x = ks[ndim*n1];
		k1.y = ks[ndim*n1+1];
		if(ndim>2)
			k1.z = ks[ndim*n1+2];
		else
			k1.z = 0.0f;
					
		for(i2=tx; i2<BINSIZE128; i2++)
		{
			n2 = bins[BINSIZE128*b2+i2];
		
			if(n2>=0)
			{
				k2.x = ks[ndim*n2];
				k2.y = ks[ndim*n2+1];
				if(ndim>2)
					k2.z = ks[ndim*n2+2];
				else
					k2.z = 0.0f;
					
				if(n1==n2)
				{
// 						if(Wout[n1]>0)
// 						if(n1==0)
// 							printf("Wout %.0f, b1 %d, b2 %d, tx %d, n1 %d, n2 %d, index %d\n", Wout[n1], b1, b2, tx, n1, n2, index);
					atomicAdd(&Wout[n1], Win[n1]);
// 						atomicAdd(&(Wout[n1]), 1);
// 						Wout[n1] = Win[n1];
				}
				else
				{
// 						if(n1==60732 || n2==60732)
// 							printf("Wout %.0f, b1 %d, b2 %d, tx %d, n1 %d, n2 %d, index %d, k1 (%f,%f,%f) k2 (%f,%f,%f)\n", Wout[n1], b1, b2, tx, n1, n2, index, k1.x, k1.y, k1.z, k2.x, k2.y, k2.z);
					
					d2 = (k2.x-k1.x)*(k2.x-k1.x) + (k2.y-k1.y)*(k2.y-k1.y) + (k2.z-k1.z)*(k2.z-k1.z);
					if(d2<=(kernelRad*kernelRad))
					{
						d = sqrtf(d2);
// 							P = P2D(d, F);
						P = P3D(d, F);
// 							P = 0.0f;
						atomicAdd(&Wout[n1], P*Win[n2]);
						atomicAdd(&Wout[n2], P*Win[n1]);
						/*atomicAdd(&(Wout[n1]), 1);
						atomicAdd(&(Wout[n2]), 1);*/
					}
				}
			}
			else
				break;
		}
	}
	else
	{
		k1.x = ks[ndim*n1];
		k1.y = ks[ndim*n1+1];
		if(ndim>2)
			k1.z = ks[ndim*n1+2];
		else
			k1.z = 0.0f;
	
		for(i2=0; i2<BINSIZE128; i2++)
		{
			n2 = bins[BINSIZE128*b2+i2];
		
			if(n2>=0)
			{
				k2.x = ks[ndim*n2];
				k2.y = ks[ndim*n2+1];
				if(ndim>2)
					k2.z = ks[ndim*n2+2];
				else
					k2.z = 0.0f;
						
				d2 = (k2.x-k1.x)*(k2.x-k1.x) + (k2.y-k1.y)*(k2.y-k1.y) + (k2.z-k1.z)*(k2.z-k1.z);
				if(d2<=(kernelRad*kernelRad))
				{
					d = sqrtf(d2);
// 						P = P2D(d, F);
					P = P3D(d, F);
// 						P = 0.0f;
					atomicAdd(&(Wout[n1]), P*Win[n2]);
					atomicAdd(&(Wout[n2]), P*Win[n1]);
					/*atomicAdd(&(Wout[n1]), 1);
					atomicAdd(&(Wout[n2]), 1);*/
// 						if(P>1.0f)
// 						if(kernelRad-1e-4<d)
// 							printf("Wout %.2f, P %f, b1 %d, b2 %d, tx %d, n1 %d, n2 %d, index %d\n", Wout[n1], P, b1, b2, tx, n1, n2, index);
				}
			}
			else
				break;
		}
	}
	
	return;
}

__global__ void conv_kernel2(unsigned int *connections, int numConnections, int *bins, float *ks, int ndim, float kernelRad, float F, float* Win, float *Wout)
{
	unsigned int index = blockIdx.x + gridDim.x*blockIdx.y;
	unsigned int tx = threadIdx.x;
	unsigned int b1, b2;
	float d, d2;
	float P;
	int i2;
	
	if(index>=numConnections)
		return;
	
	b1 = connections[2*index];
	b2 = connections[2*index+1];
	
	float3 k1;
	int n1;
	
	__shared__ float3 k2[BINSIZE128];
	__shared__ int n2[BINSIZE128];
	
	n1 = bins[BINSIZE128*b1+tx];
	
	if(n1<0)
		return;
	
	k1.x = ks[ndim*n1];
	k1.y = ks[ndim*n1+1];
	if(ndim>2)
		k1.z = ks[ndim*n1+2];
	else
		k1.z = 0.0f;
	
	n2[tx] = bins[BINSIZE128*b2+tx];
	if(n2>=0)
	{
		k2[tx].x = ks[ndim*n2[tx]];
		k2[tx].y = ks[ndim*n2[tx]+1];
		if(ndim>2)
			k2[tx].z = ks[ndim*n2[tx]+2];
		else
			k2[tx].z = 0.0f;
	}
	
	__syncthreads();
	
	if(b1==b2)
	{
		for(i2=tx; i2<BINSIZE128; i2++)
		{
			if(n2[i2]>=0)
			{
				if(n1==n2[i2])
				{
// 						if(Wout[n1]>0)
// 						if(n1==0)
// 							printf("Wout %.0f, b1 %d, b2 %d, tx %d, n1 %d, n2 %d, index %d\n", Wout[n1], b1, b2, tx, n1, n2, index);
					atomicAdd(&Wout[n1], Win[n1]);
// 						atomicAdd(&(Wout[n1]), 1);
// 						Wout[n1] = Win[n1];
				}
				else
				{
// 						if(n1==60732 || n2==60732)
// 							printf("Wout %.0f, b1 %d, b2 %d, tx %d, n1 %d, n2 %d, index %d, k1 (%f,%f,%f) k2 (%f,%f,%f)\n", Wout[n1], b1, b2, tx, n1, n2, index, k1.x, k1.y, k1.z, k2.x, k2.y, k2.z);
					
					d2 = (k2[i2].x-k1.x)*(k2[i2].x-k1.x) + (k2[i2].y-k1.y)*(k2[i2].y-k1.y) + (k2[i2].z-k1.z)*(k2[i2].z-k1.z);
					if(d2<=(kernelRad*kernelRad))
					{
						d = sqrtf(d2);
						P = P2D(d, F);
						atomicAdd(&Wout[n1], P*Win[n2[i2]]);
						atomicAdd(&Wout[n2[i2]], P*Win[n1]);
						/*atomicAdd(&(Wout[n1]), 1);
						atomicAdd(&(Wout[n2]), 1);*/
					}
				}
			}
			else
				break;
		}
	}
	else
	{
		for(i2=0; i2<BINSIZE128; i2++)
		{
			if(n2[i2]>=0)
			{
				d2 = (k2[i2].x-k1.x)*(k2[i2].x-k1.x) + (k2[i2].y-k1.y)*(k2[i2].y-k1.y) + (k2[i2].z-k1.z)*(k2[i2].z-k1.z);
				if(d2<=(kernelRad*kernelRad))
				{
					d = sqrtf(d2);
					P = P2D(d, F);
					atomicAdd(&(Wout[n1]), P*Win[n2[i2]]);
					atomicAdd(&(Wout[n2[i2]]), P*Win[n1]);
					/*atomicAdd(&(Wout[n1]), 1);
					atomicAdd(&(Wout[n2]), 1);*/
// 						if(P>1.0f)
// 						if(kernelRad-1e-4<d)
// 							printf("Wout %.2f, P %f, b1 %d, b2 %d, tx %d, n1 %d, n2 %d, index %d\n", Wout[n1], P, b1, b2, tx, n1, n2, index);
				}
			}
			else
				break;
		}
	}
	
	return;
}

/**
Using template
*/
template <int NDIM, int BINSIZE>
__global__ void conv_kernelt(unsigned int *connections, int numConnections, int *bins, float *ks, float kernelRad, float F, float* Win, float *Wout)
{
	unsigned int index = blockIdx.x + gridDim.x*blockIdx.y;
	unsigned int tx = threadIdx.x;
	unsigned int b1, b2;
	float3 k1, k2;
	float d, d2;
	float P;
	int n1, n2;
	int i2;
	
	if(index>=numConnections)
		return;
	
	b1 = connections[2*index];
	b2 = connections[2*index+1];
	
	n1 = bins[BINSIZE*b1+tx];
	if(n1<0)
		return;
	
	if(b1==b2)
	{
//		Processing matches in the same bin
		k1.x = ks[NDIM*n1];
		k1.y = ks[NDIM*n1+1];
		if(NDIM>2)
			k1.z = ks[NDIM*n1+2];
		else
			k1.z = 0.0f;
					
		for(i2=tx; i2<BINSIZE; i2++)
		{
			n2 = bins[BINSIZE*b2+i2];
		
			if(n2>=0)
			{
				k2.x = ks[NDIM*n2];
				k2.y = ks[NDIM*n2+1];
				if(NDIM>2)
					k2.z = ks[NDIM*n2+2];
				else
					k2.z = 0.0f;
					
				if(n1==n2)
				{
// 						if(Wout[n1]>0)
// 						if(n1==0)
// 							printf("Wout %.0f, b1 %d, b2 %d, tx %d, n1 %d, n2 %d, index %d\n", Wout[n1], b1, b2, tx, n1, n2, index);
					atomicAdd(&Wout[n1], Win[n1]);
// 						atomicAdd(&(Wout[n1]), 1);
// 						Wout[n1] = Win[n1];
				}
				else
				{
// 						if(n1==60732 || n2==60732)
// 							printf("Wout %.0f, b1 %d, b2 %d, tx %d, n1 %d, n2 %d, index %d, k1 (%f,%f,%f) k2 (%f,%f,%f)\n", Wout[n1], b1, b2, tx, n1, n2, index, k1.x, k1.y, k1.z, k2.x, k2.y, k2.z);
					
					if(NDIM==3)
						d2 = (k2.x-k1.x)*(k2.x-k1.x) + (k2.y-k1.y)*(k2.y-k1.y) + (k2.z-k1.z)*(k2.z-k1.z);
					else
						d2 = (k2.x-k1.x)*(k2.x-k1.x) + (k2.y-k1.y)*(k2.y-k1.y);

					if(d2<=(kernelRad*kernelRad))
					{
						d = sqrtf(d2);
						if(NDIM==3)
							P = P3D(d, F);
						else
 							P = P2D(d, F);
// 							P = 0.0f;
						atomicAdd(&Wout[n1], P*Win[n2]);
						atomicAdd(&Wout[n2], P*Win[n1]);
						/*atomicAdd(&(Wout[n1]), 1);
						atomicAdd(&(Wout[n2]), 1);*/
					}
				}
			}
			else
				break;
		}
	}
	else
	{
		k1.x = ks[NDIM*n1];
		k1.y = ks[NDIM*n1+1];
		if(NDIM>2)
			k1.z = ks[NDIM*n1+2];
		else
			k1.z = 0.0f;
	
		for(i2=0; i2<BINSIZE; i2++)
		{
			n2 = bins[BINSIZE*b2+i2];
		
			if(n2>=0)
			{
				k2.x = ks[NDIM*n2];
				k2.y = ks[NDIM*n2+1];
				if(NDIM>2)
					k2.z = ks[NDIM*n2+2];
				else
					k2.z = 0.0f;
						
				d2 = (k2.x-k1.x)*(k2.x-k1.x) + (k2.y-k1.y)*(k2.y-k1.y); 
				if(NDIM==3)
					d2 += (k2.z-k1.z)*(k2.z-k1.z);

				if(d2<=(kernelRad*kernelRad))
				{
					d = sqrtf(d2);
					if(NDIM==2)
 						P = P2D(d, F);
					else
						P = P3D(d, F);
// 						P = 0.0f;
					atomicAdd(&(Wout[n1]), P*Win[n2]);
					atomicAdd(&(Wout[n2]), P*Win[n1]);
					/*atomicAdd(&(Wout[n1]), 1);
					atomicAdd(&(Wout[n2]), 1);*/
// 						if(P>1.0f)
// 						if(kernelRad-1e-4<d)
// 							printf("Wout %.2f, P %f, b1 %d, b2 %d, tx %d, n1 %d, n2 %d, index %d\n", Wout[n1], P, b1, b2, tx, n1, n2, index);
				}
			}
			else
				break;
		}
	}

	//atomicAdd(pct, 1.0f/numConnections);
	//printf("%.2f\t", *pct);
	
	return;
}

/**
Separable kernel
*/
template <int NDIM, int BINSIZE>
__global__ void conv_kernel_sep(unsigned int *connections, int numConnections, int *bins, float *ks, float* kernelRad, float *F, float* Win, float *Wout)
{
	unsigned int index = blockIdx.x + gridDim.x*blockIdx.y;
	unsigned int tx = threadIdx.x;
	unsigned int b1, b2;
	float3 k1, k2;
	float3 d;
	float P;
	int n1, n2;
	int i2;
	
	if(index>=numConnections)
		return;
	
	b1 = connections[2*index];
	b2 = connections[2*index+1];
	
	n1 = bins[BINSIZE*b1+tx];
	if(n1<0)
		return;
		
	if(b1==b2)
	{
		k1.x = ks[NDIM*n1];
		k1.y = ks[NDIM*n1+1];
		if(NDIM>2)
			k1.z = ks[NDIM*n1+2];
		else
			k1.z = 0.0f;
					
		for(i2=tx; i2<BINSIZE; i2++)
		{
			n2 = bins[BINSIZE*b2+i2];
		
			if(n2>=0)
			{
				k2.x = ks[NDIM*n2];
				k2.y = ks[NDIM*n2+1];
				if(NDIM>2)
					k2.z = ks[NDIM*n2+2];
				else
					k2.z = 0.0f;
					
				if(n1==n2)
				{
// 						if(Wout[n1]>0)
// 						if(n1==0)
// 							printf("Wout %.0f, b1 %d, b2 %d, tx %d, n1 %d, n2 %d, index %d\n", Wout[n1], b1, b2, tx, n1, n2, index);
					atomicAdd(&Wout[n1], Win[n1]);
// 						atomicAdd(&(Wout[n1]), 1);
// 						Wout[n1] = Win[n1];
				}
				else
				{
// 						if(n1==60732 || n2==60732)
// 							printf("Wout %.0f, b1 %d, b2 %d, tx %d, n1 %d, n2 %d, index %d, k1 (%f,%f,%f) k2 (%f,%f,%f)\n", Wout[n1], b1, b2, tx, n1, n2, index, k1.x, k1.y, k1.z, k2.x, k2.y, k2.z);
					d.x = fabs(k2.x-k1.x);
					d.y = fabs(k2.y-k1.y);
					d.z = fabs(k2.z-k1.z);

					if(NDIM==2)
					{
						if((d.x<kernelRad[0]) && (d.y<kernelRad[1]))
						{
							P = P1D(d.x, F[0]) * P1D(d.y, F[1]);
							atomicAdd(&Wout[n1], P*Win[n2]);
							atomicAdd(&Wout[n2], P*Win[n1]);
						}
					}
					else
					{
						if((d.x<kernelRad[0]) && (d.y<kernelRad[1]) && (d.z<kernelRad[2]))
						{
							P = P1D(d.x, F[0]) * P1D(d.y, F[1]) * P1D(d.z, F[2]);
							atomicAdd(&Wout[n1], P*Win[n2]);
							atomicAdd(&Wout[n2], P*Win[n1]);
						}	
					}
					
						/*atomicAdd(&(Wout[n1]), 1);
						atomicAdd(&(Wout[n2]), 1);*/
				}
			}
			else
				break;
		}
	}
	else
	{
		k1.x = ks[NDIM*n1];
		k1.y = ks[NDIM*n1+1];
		if(NDIM>2)
			k1.z = ks[NDIM*n1+2];
		else
			k1.z = 0.0f;
	
		for(i2=0; i2<BINSIZE; i2++)
		{
			n2 = bins[BINSIZE*b2+i2];
		
			if(n2>=0)
			{
				k2.x = ks[NDIM*n2];
				k2.y = ks[NDIM*n2+1];
				if(NDIM>2)
					k2.z = ks[NDIM*n2+2];
				else
					k2.z = 0.0f;
						
				
					d.x = fabs(k2.x-k1.x);
					d.y = fabs(k2.y-k1.y);
					d.z = fabs(k2.z-k1.z);

					if(NDIM==2)
					{
						if((d.x<kernelRad[0]) && (d.y<kernelRad[1]))
						{
							P = P1D(d.x, F[0]) * P1D(d.y, F[1]);
							atomicAdd(&Wout[n1], P*Win[n2]);
							atomicAdd(&Wout[n2], P*Win[n1]);
						}
					}
					else
					{
						if((d.x<kernelRad[0]) && (d.y<kernelRad[1]) && (d.z<kernelRad[2]))
						{
							P = P1D(d.x, F[0]) * P1D(d.y, F[1]) * P1D(d.z, F[2]);
							atomicAdd(&Wout[n1], P*Win[n2]);
							atomicAdd(&Wout[n2], P*Win[n1]);
						}	
					}

					/*atomicAdd(&(Wout[n1]), 1);
					atomicAdd(&(Wout[n2]), 1);*/
// 						if(P>1.0f)
// 						if(kernelRad-1e-4<d)
// 							printf("Wout %.2f, P %f, b1 %d, b2 %d, tx %d, n1 %d, n2 %d, index %d\n", Wout[n1], P, b1, b2, tx, n1, n2, index);
			}
			else
				break;
		}
	}
	
	return;
}

__global__ void divideArray(float* num, float* den, int npts)
{
	unsigned int index = threadIdx.x + blockDim.x*blockIdx.x;
	
	if(index>=npts)
		return;
	
	if(fabsf(den[index])>0.0f)
		num[index] /= den[index];
	else
		num[index] = 0.0f;
	
	return;
}

int findDataSize(int target, int max1, int max2, unsigned int *dim1, unsigned int *dim2)
{
     int success= 0;
     int rem;
     int minRem=0;
     int d1, d2;

     if(target<max1)
     {
          *dim1 = target;
          *dim2 = 1;
     }
     else
     {
          d1 = max1;

          while(d1>0)
          {
               d2 = ceil(target/(1.0*d1));
               rem = d1*d2 - target;

               if(d2>max2)
                    d1 = 0;
               else if(rem<=minRem)
               {
                    *dim1 = d1;
                    *dim2 = d2;
                    success = 1;
               }
               d1--;
          }
     }

     return success;
}

/**
 \param[in]	ndim	Number of dimensions
 \param[in]	FOV	Field of View (cm)
*/
void jdcfcu(int npts, float *k, int ndim, int *nInclude, float *FOV, float *voxDim, ConvType convType, int numLobes, int *numCt, int compartmentSize, int numIter, float *dcf)
{
	unsigned int *binIdx;
	unsigned int *binCount1;
	int *bins, *bins_dev;
	unsigned int *binStartIdx;
	unsigned int *binCurrentIdx;
	unsigned int *binReps;
	unsigned int *binConnections, *binConnections_dev;
	int numConnections;
	int numBins1 = 1;
	int numBins2 = 0;
	int bsi, bci;
	int bi, bii;
	int freeBinIdx;
	int c[3];
	int nptsActual;
	float pctBinFull;
	float kmax[3];
	float F[3];
	float *F_dev;
	float kernelRad[3];
	float *kernelRad_dev;
	float kernelRadMax;
	int maxCtOffset[3] = {1,1,1};
	float *Win_dev;
	float *Wout_dev;
	float *k_dev;
	dim3 grid_size, block_size;
	dim3 grid_size_div, block_size_div;
	
	int dim;
	int m, n, p, q;
	int mm, nn, pp, qq;
	int nnStart, ppStart, qqStart;
	
	printf("\n%d iterations\n", numIter);
	printf("%d lobes\n", numLobes);
	printf("Compartment size %d\n", compartmentSize);

	//Setup kernels
	for(n=0; n<ndim; n++)
	{
		kmax[n] = 5/voxDim[n];
		F[n] = FOV[n]*2*kmax[n];	/* parameter for kernel width */
	}

	switch(convType)
	{
		case ctCONST:
			if(ndim==2)
				*kernelRad = (numLobes*1.01+0.21)/(*F);
			else
				*kernelRad = (numLobes*1.01+0.44f)/(*F);
			kernelRadMax = *kernelRad;

			break;
		case ctSEP:
			kernelRadMax = 0.0f;
			for(n=0; n<ndim; n++)
			{
				kernelRad[n] = numLobes/F[n];	// sinc has spacing 1
				kernelRadMax = max(kernelRadMax, kernelRad[n]);
			}

			break;
	}

//	Calculate optimal bin size
	for(dim=0; dim<ndim; dim++)
	{
		if(numCt[dim]<0)
		{
			switch(convType)
			{
				case ctCONST:
					numCt[dim] = 1/kernelRad[0]-1;
					break;
				case ctSEP:
					numCt[dim] = 1/kernelRad[dim]-1;
					break;
			}
			numCt[dim] -= numCt[dim]%2;
		}
	}

//	Calculate total number of bins
	for(dim=0; dim<ndim; dim++)
		numBins1 *= numCt[dim];
	
	binCount1 = (unsigned int*)malloc(numBins1*sizeof(unsigned int));

	// Initialize number of points in grid to 0
// 	memste(binCount1, 0, numBins1*sizeof(unsigned int));
	for(n=0; n<numBins1; n++)
		binCount1[n] = 0;
		
	for(dim=ndim; dim<3; dim++)
		numCt[dim] = 1;
	
	if(nInclude==NULL)
	{
		nInclude = (int*)malloc(npts*sizeof(int));
		for(n=0; n<npts; n++)
			nInclude[n] = 1;
	}
	
	binIdx = (unsigned int*)malloc(npts*sizeof(unsigned int));
	nptsActual = 0;
	for(n=0; n<npts; n++)
	{
		if(nInclude[n])
		{
// 			dcf[n] = 1.0f;
			for(dim=0; dim<ndim; dim++)
				c[dim] = floor(k[ndim*n+dim]*numCt[dim]+numCt[dim]/2.0);
			for(dim=ndim; dim<3; dim++)
				c[dim] = 0;
			
			/* Make sure c is in bounds */
			for(dim=0; dim<ndim; dim++)
			{
				c[dim] = min(c[dim], numCt[dim]-1);
				c[dim] = max(c[dim], 0);
			}
				
			binIdx[n] = (c[2]*numCt[1]+c[1])*numCt[0] + c[0];
			if(binIdx[n]>numBins1)
			{
				fprintf(stderr, "Error numbins=%d, current bin idx=%d, n= %d\n", numBins1, binIdx[n], n);
				return;
			}
			binCount1[binIdx[n]]++;
			nptsActual++;
		}
		else
			dcf[n] = 0.0f;
	}
	
	printf("Calculating for %d/%d points\n", nptsActual, npts);

	//Count number of required bins
	binStartIdx = (unsigned int*)malloc(numBins1*sizeof(unsigned int));
	binReps = (unsigned int*)malloc(numBins1*sizeof(unsigned int));
	for(n=0; n<numBins1; n++)
	{
		if(n==0)
			binStartIdx[n] = 0;
		else
			binStartIdx[n] = binStartIdx[n-1]+binReps[n-1];
		
		binReps[n] = ceil(binCount1[n]/(float)compartmentSize);
		numBins2 += binReps[n];
	}

	bins = (int*)malloc(compartmentSize *numBins2*sizeof(int));

	binCurrentIdx = (unsigned int*)malloc(numBins1*sizeof(unsigned int));
	memset(binCurrentIdx, 0, numBins1*sizeof(unsigned int));

	for(n=0; n< compartmentSize *numBins2; n++)
		bins[n] = -1;

	for(n=0; n<npts; n++)
	{
		if(nInclude[n])
		{
			bi = binIdx[n];
			bci = binCurrentIdx[bi];
			bsi = compartmentSize *binStartIdx[bi];
			bins[bsi+bci++] = n;
			binCurrentIdx[bi] = bci;
		}
	}

	pctBinFull = 0;
	for(n=0; n< compartmentSize *numBins2; n++)
		if(bins[n]!=-1)
			pctBinFull += 1.0f/(compartmentSize *numBins2);

	printf("%d bins\n", numBins2);
	printf("%f%% of binspace is full\n", 100*pctBinFull);
	
// 	Setup connections
	if(convType==ctCONST)
	{
		printf("Kernel radius %f\n", *kernelRad);
		printf("F %f\n", *F);
	}
	else
	{
		printf("Kernel radii %f %f %f\n", kernelRad[0], kernelRad[1], kernelRad[2]);
		printf("F %f %f %f\n", F[0], F[1], F[2]);
	}

//	Calculate maximum compartment search distance
	switch(convType)
	{
		case ctCONST:
			for(dim=0; dim<ndim; dim++)
				maxCtOffset[dim] = (int)(*kernelRad*numCt[dim]+1);
			break;
		case ctSEP:
			for(n=0; n<ndim; n++)
				maxCtOffset[n] = (int)(kernelRad[n]*numCt[n]+1);
			break;
	}

	
/*	if(ndim==2)
		binConnections = (unsigned int*)malloc(2*100000*sizeof(unsigned int));
	else if(ndim==3)
		binConnections = (unsigned int*)malloc(2*20000000*sizeof(unsigned int));*/
		
	for(dim=0; dim<2; dim++)
	{
		// on first pass count number of connections and allocate
		// on second pass assign connections
		if(dim==1)
			binConnections = (unsigned int*)malloc(2*numConnections*sizeof(unsigned int));
		numConnections = 0;

	for(m=0; m<numCt[0]; m++)
		for(n=0; n<numCt[1]; n++)
			for(p=0; p<numCt[2]; p++)
			{
				bi = (p*numCt[1]+n)*numCt[0]+m;
				
				for(mm=m; mm<min(m+maxCtOffset[0]+1, numCt[0]); mm++)
				{
					if(mm==m)
						nnStart = n;
					else
						nnStart = max(0,n-maxCtOffset[1]);
					for(nn=nnStart; nn<min(n+maxCtOffset[1]+1, numCt[1]); nn++)
					{
						if((nn==n) && (mm==m))
							ppStart = p;
						else
							ppStart = max(0, p-maxCtOffset[2]);
									
						for(pp=ppStart; pp<min(p+maxCtOffset[2]+1, numCt[2]); pp++)
						{
							bii = (pp*numCt[1]+nn)*numCt[0]+mm;
							for(q=0; q<binReps[bi]; q++)
							{
								if(bi==bii)
									qqStart = q;
								else
									qqStart = 0;
								for(qq=qqStart; qq<binReps[bii]; qq++)
								{
									if(dim==0)
										numConnections++;
									else
									{
										binConnections[2*numConnections] = binStartIdx[bi]+q;
										binConnections[2*numConnections+++1] = binStartIdx[bii]+qq;
									}
// 									printf("numConnections %d\n", numConnections);
								}
							}
						}
					}
				}
			}
	}
	
	printf("# compartments  %d %d %d\n", numCt[0], numCt[1], numCt[2]);
	printf("Max compartment offset  %d %d %d\n", maxCtOffset[0], maxCtOffset[1], maxCtOffset[2]);
	printf("%d connections\n", numConnections);
	
	block_size.x = compartmentSize;
// 	block_size.y = BINSIZE;
	grid_size.x = numConnections;
	if(!findDataSize(numConnections, 65536, 65536, &grid_size.x, &grid_size.y))
	{
// 		Coldn't find exact factors
		grid_size.x = sqrt(numConnections);
		grid_size.y =ceil(numConnections/(double)grid_size.x);
	}
	printf("Grid size [%d,%d,%d]\n", grid_size.x, grid_size.y, grid_size.z);
	printf("Block size [%d,%d,%d]\n", block_size.x, block_size.y, block_size.z);
	
 	findBlockGrid3(npts, 256, &block_size_div, &grid_size_div);
	block_size_div.x = 256;
/*	grid_size_div.x = npts/block_size_div.x;
	grid_size_div.x++;*/
	printf("Divide grid size [%d,%d,%d]\n", grid_size_div.x, grid_size_div.y, grid_size_div.z);
	printf("Divide block size [%d,%d,%d]\n", block_size_div.x, block_size_div.y, block_size_div.z);
	
// 	Allocate gpu memory
  timer.start();
	cudaMalloc((void**)&binConnections_dev, 2*numConnections*sizeof(unsigned int));
	cudaMalloc((void**)&bins_dev, compartmentSize*numBins2*sizeof(int));
	cudaMalloc((void**)&k_dev, ndim*npts*sizeof(float));
	cudaMalloc((void**)&Win_dev, npts*sizeof(float));
	cudaMalloc((void**)&Wout_dev, npts*sizeof(float));
	
// 	Copy to gpu
	cudaMemcpy(binConnections_dev, binConnections, 2*numConnections*sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(bins_dev, bins, compartmentSize*numBins2*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(k_dev, k, ndim*npts*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(Win_dev, dcf, npts*sizeof(float), cudaMemcpyHostToDevice);
  timer.stop("Allocate and copy to gpu");
	
	for(n=0; n<numIter; n++)
	{
    timer.start();
		cudaMemset(Wout_dev, 0, npts*sizeof(float));
//		conv_kernel<<<grid_size, block_size>>>(binConnections_dev, numConnections, bins_dev, k_dev, ndim, kernelRad, F, Win_dev, Wout_dev);
		switch(convType)
		{
			case ctCONST:
				if(ndim==2)
					conv_kernelt<2, 128><<<grid_size, block_size>>>(binConnections_dev, numConnections, bins_dev, k_dev, *kernelRad, *F, Win_dev, Wout_dev);
				else
          switch (compartmentSize)
          {
            case COMPARTMENT_128:
					    conv_kernelt<3, COMPARTMENT_128><<<grid_size, block_size>>>(binConnections_dev, numConnections, bins_dev, k_dev, *kernelRad, *F, Win_dev, Wout_dev);
              break;
            case COMPARTMENT_64:
              conv_kernelt<3, COMPARTMENT_64><<<grid_size, block_size>>>(binConnections_dev, numConnections, bins_dev, k_dev, *kernelRad, *F, Win_dev, Wout_dev);
              break;
          }
				break;
			case ctSEP:
				cudaMalloc((void**)&kernelRad_dev, 3*sizeof(float));
				cudaMalloc((void**)&F_dev, 3*sizeof(float));

				cudaMemcpy(kernelRad_dev, kernelRad, 3*sizeof(float), cudaMemcpyHostToDevice);
				cudaMemcpy(F_dev, F, 3*sizeof(float), cudaMemcpyHostToDevice);
				checkLaunch("copy");

				if(ndim==2)
					switch(compartmentSize)
					{
						case COMPARTMENT_64:
							conv_kernel_sep<2, COMPARTMENT_64><<<grid_size, block_size>>>(binConnections_dev, numConnections, bins_dev, k_dev, kernelRad_dev, F_dev, Win_dev, Wout_dev);
							break;
						case COMPARTMENT_128:
							conv_kernel_sep<2, COMPARTMENT_128><<<grid_size, block_size>>>(binConnections_dev, numConnections, bins_dev, k_dev, kernelRad_dev, F_dev, Win_dev, Wout_dev);
							break;
					}
				else
					switch(compartmentSize)
					{
						case COMPARTMENT_64:
							conv_kernel_sep<3, COMPARTMENT_64><<<grid_size, block_size>>>(binConnections_dev, numConnections, bins_dev, k_dev, kernelRad_dev, F_dev, Win_dev, Wout_dev);
							break;
						case COMPARTMENT_128:
							conv_kernel_sep<3, COMPARTMENT_128><<<grid_size, block_size>>>(binConnections_dev, numConnections, bins_dev, k_dev, kernelRad_dev, F_dev, Win_dev, Wout_dev);
							break;
						case 256:
							conv_kernel_sep<3, 256><<<grid_size, block_size>>>(binConnections_dev, numConnections, bins_dev, k_dev, kernelRad_dev, F_dev, Win_dev, Wout_dev);
							break;
						case 512:
							conv_kernel_sep<3, 512><<<grid_size, block_size>>>(binConnections_dev, numConnections, bins_dev, k_dev, kernelRad_dev, F_dev, Win_dev, Wout_dev);
							break;
						default:
							fprintf(stderr, "Unsupported bin size %d {64,128,256,512}\n", compartmentSize);
							abort();
					}
				break;
		}
		checkLaunch("convolution");
		divideArray<<<grid_size_div, block_size_div>>>(Win_dev, Wout_dev, npts);
		checkLaunch("division");
    timer.stop("convolution");
	}
	
	
	cudaMemcpy(dcf, Wout_dev, npts*sizeof(float), cudaMemcpyDeviceToHost);
	sprintf(filename, "%sWout", savePath);
	file = fopen(filename, "wb");
	if(file!=NULL)
	{
		fwrite(dcf, npts, sizeof(float), file);
		fclose(file);
	}
	
	cudaMemcpy(dcf, Win_dev, npts*sizeof(float), cudaMemcpyDeviceToHost);
	
	sprintf(filename, "%sbinConnections", savePath);
	file = fopen(filename, "wb");
	if(file!=NULL)
	{
		fwrite(binConnections, 2*numConnections, sizeof(unsigned int), file);
		fclose(file);
	}
	
	sprintf(filename, "%sbins", savePath);
	file = fopen(filename, "wb");
	if(file!=NULL)
	{
		fwrite(bins, compartmentSize*numBins2, sizeof(int), file);
		fclose(file);
	}
	
	cudaFree(binConnections_dev);
	cudaFree(bins_dev);
	cudaFree(Wout_dev);
	cudaFree(Win_dev);
	cudaFree(k_dev);
	
	free(binIdx);
	free(binConnections);

	return;
}

#ifndef MATLAB_MEX_FILE

extern "C"{
//#include "cones.h"
//#include "spiral.h"
}

int main(int argc, char *argv[])
{
	int numLobes = 2;
	int numCt[] = {-1,-1,-1};
	int numIter = 1;
	ConvType cType = ctCONST;

	//printTrajectoryInfo(&traj, NULL);

  const char* filePathInput = argv[1];
  printf("Loading k-space file from %s\n", filePathInput);

  Endian endian = argv[11][0] == 'b' ? BigEndian : LittleEndian;
  Trajectory* trajectory = loadKSpaceFile(filePathInput, atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), endian);
  std::vector<float> coordinates;
  std::vector<float> densityCompensation;
	for(int r=0; r<trajectory->numReadouts; r++)
			for(int n=0; n< trajectory->numReadoutPoints; n++)
			{
        float pointCoordinates[3];
        float density;
        trajectoryCoordinates(n, r, trajectory, pointCoordinates, &density);
        for (int d = 0; d < trajectory->numDimensions; d++)
          coordinates.push_back(pointCoordinates[d]);
        densityCompensation.push_back(density);
			}
	
  float fieldOfView[3];
  float spatialResolution[3];
  const int fieldOfViewArgumentIndex = 5;
  for (int d = 0; d < trajectory->numDimensions; d++)
  {
    fieldOfView[d] = atoi(argv[d + fieldOfViewArgumentIndex]);
    spatialResolution[d] = atoi(argv[d + fieldOfViewArgumentIndex + trajectory->numDimensions]);
  }
  jdcfcu(densityCompensation.size(), coordinates.data(), trajectory->numDimensions, NULL, fieldOfView, spatialResolution, cType, numLobes, numCt, COMPARTMENT_64, numIter, densityCompensation.data());
  for(int n=0; n<densityCompensation.size(); n++)
  {
    const int readout = n / trajectory->numReadoutPoints;
    const int point = n % trajectory->numReadoutPoints;
    setTrajectoryPoint(point, readout, trajectory, NULL, densityCompensation[n]);
  }

  const char* filePathOutput = argv[12];
  printf("Writing new file %s\n", filePathOutput);
  saveKSPaceFile(filePathOutput, trajectory, endian);
}

#else

#include "mex.h"

/*
0: k-space trajectory [axis, readout, interleaf]
1: dcf mask
2: FOV (cm)
3: Voxel dimensions (mm)
4: Convolution type 0-Constant 1-Separable
5: # of lobes
6: # of compartments
7: # iterations
8: dcf estimate
*/
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	float FOV[3];
	float vd[3];
	const mwSize *trajDims;
	int numDim;	// number of dimensions
	const mwSize *dcfDims;
	float *k;
	float *dcf;
	int *dcfInclude;	// indeces to include in processing
	int numCt[] = {-1,-1,-1};
	int numLobes = 2;
	int numIter = 1;	// number of iterations
	ConvType cType = ctSEP;
	int npts;	// number of point in trajectory
	int n;
	int doLoadDcfMask = 1;
	int doLoadDcf = 0;
	
	double *dptr;
	char msg[128];
	
	if(!mxIsSingle(prhs[0]))
		mexErrMsgTxt("k-space trajectory not given as single precission float");
	
	trajDims = mxGetDimensions(prhs[0]);
	numDim = trajDims[0];
	if(numDim>3)
		mexErrMsgTxt("More than 3 axes given in trajectory");
	
	npts = trajDims[1]*trajDims[2];
	
	k = (float*)mxGetData(prhs[0]);
	
// 	dcf mask
	if(nrhs<2)
		mexErrMsgTxt("dcf mask not given as 4 byte int");
	
	if(mxIsEmpty(prhs[1]))
	{
// 		if dcf mask not specified include every point
		doLoadDcfMask = 0;
		dcfInclude = (int*)malloc(npts*sizeof(int));
		for(n=0; n<npts; n++)
			dcfInclude[n] = 1;
	}
	else
	{
		if(!mxIsInt32(prhs[1]))
			mexErrMsgTxt("dcf mask not given as 4 byte int");
		dcfInclude = (int*)mxGetData(prhs[1]);
	}
	
// 	FOV
	if(nrhs<3)
		mexErrMsgTxt("FOV not given");
	else if(mxGetNumberOfElements(prhs[2])<numDim)
	{
		sprintf(msg, "%d dimension(s) of FOV given, need %d", mxGetNumberOfElements(prhs[2]), numDim);
		mexErrMsgTxt(msg);
	}
	else
	{
		dptr = mxGetPr(prhs[2]);
		for(n=0; n<trajDims[0]; n++)
			FOV[n] = dptr[n];
	}
	
// 	Resolution
	if(nrhs<4)
		mexErrMsgTxt("resolution not given");
	else if(mxGetNumberOfElements(prhs[3])<numDim)
	{
		sprintf(msg, "%d voxel dimension(s) given, need %d", mxGetNumberOfElements(prhs[3]), numDim);
		mexErrMsgTxt(msg);
	}
	else
	{
		dptr = mxGetPr(prhs[3]);
		for(n=0; n<trajDims[0]; n++)
			vd[n] = dptr[n];
	}
	
// 	Convolution type
	if(nrhs>4)
		cType = (ConvType)mxGetScalar(prhs[4]);

// 	Number of lobes
	if(nrhs>5)
		numLobes = mxGetScalar(prhs[5]);

	if(nrhs>6)
	{
		if(mxGetNumberOfElements(prhs[6])<numDim)
		{
			sprintf(msg, "%d number compartment(s) given, need %d", mxGetNumberOfElements(prhs[6]), numDim);
			mexErrMsgTxt(msg);
		}
		dptr = mxGetPr(prhs[6]);
		for(n=0; n<trajDims[0]; n++)
			numCt[n] = dptr[n];
	}
		
	if(nrhs>7)
		numIter = mxGetScalar(prhs[7]);
	
	if(DBG_DCF)
	{
		mexPrintf("%d points\n", npts);
		mexPrintf("FOV %f %f %f\n", FOV[0], FOV[1], FOV[2]);
		mexPrintf("# lobes %d\n", numLobes);
		mexPrintf("# compartments %d %d %d\n", numCt[0], numCt[1], numCt[2]);
		mexPrintf("# iterations %d\n", numIter);
	}
	
	if(nrhs>8)
	{
// 		Update the dcf given
		if(!mxIsSingle(prhs[8]))
			mexErrMsgTxt("dcf estimate not given as float");
		
		dcfDims = mxGetDimensions(prhs[8]);
		for(n=0; n<2; n++)
			if(dcfDims[n] != trajDims[n+1])
			{
				sprintf(msg, "Dcf dimensions [%d, %d] do not equal trajectory dimensions [%d, %d]\n", dcfDims[0], dcfDims[1], trajDims[1], trajDims[2]);
				mexErrMsgTxt(msg);
			}
		plhs[0] = mxDuplicateArray(prhs[8]);
		dcf = (float*)mxGetData(plhs[0]);
	}
	else
	{
// 		Start from flat dcf
		plhs[0] = mxCreateNumericArray(2, dcfDims, mxSINGLE_CLASS, mxREAL);
		dcf = (float*)mxGetData(plhs[0]);
		for(n=0; n<npts; n++)
			dcf[n] = 1.0f;
	}
	
	if(DBG_DCF) mexPrintf("Trajectory dimensions [%d %d]\n", dcfDims[0], dcfDims[1]);

	jdcfcu(npts, k, trajDims[0], dcfInclude, FOV, vd, cType, numLobes, numCt, 64, numIter, dcf);
	
	if(!doLoadDcfMask)
		free(dcfInclude);
	
	return;
}

#endif

