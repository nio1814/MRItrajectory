#include <stdio.h>
#include <string.h>

#include "mp1-util.h"

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

int checkLaunch(const std::string kernel_name)
{
	char msg[128];
	
  cudaThreadSynchronize();
  cudaError lastErr = cudaGetLastError();
  
#ifndef MATLAB_MEX_FILE
  if(lastErr == cudaSuccess)
		printf("done with %s kernel\n", kernel_name.c_str());
	else
	{
		printf("error %d '%s' on %s kernel\n", lastErr, cudaGetErrorString(lastErr), kernel_name.c_str());
		exit(1);
	}
#else
	if(lastErr == cudaSuccess)
		mexPrintf("done with %s kernel\n",kernel_name);
	else
	{
		sprintf(msg, "error %d '%s' on %s kernel\n", lastErr, cudaGetErrorString(lastErr), kernel_name);
		mexErrMsgTxt(msg);
	}
#endif
  
  return lastErr;
}

void printExit(char* msg)
{
#ifdef MATLAB_MEX_FILE
	mexErrMsgTxt(msg);
#else
	fprintf(stderr, "%s\n", msg);
	exit(EXIT_FAILURE);
#endif
}


void Timer::start()
{
	int gpuIDcurrent;	//Gpu device called from
	cudaGetDevice(&gpuIDcurrent);
	cudaSetDevice(0);
	
  cudaEventCreate(&this->startEvent);
  cudaEventCreate(&this->endEvent);
  cudaEventRecord(this->startEvent, 0);
  
	cudaSetDevice(gpuIDcurrent);
}


float Timer::stop(const std::string description)
{
	int gpuIDcurrent;	//Gpu device called from
	cudaGetDevice(&gpuIDcurrent);
	cudaSetDevice(0);
	
  cudaEventRecord(this->endEvent, 0);
  cudaEventSynchronize(this->endEvent);
  
  float elapsed_time;
  cudaEventElapsedTime(&elapsed_time, this->startEvent, this->endEvent);
  printf("%s took %.1f ms\n", description.c_str(), elapsed_time);
  cudaEventDestroy(this->startEvent);
  cudaEventDestroy(this->endEvent);
  
  cudaSetDevice(gpuIDcurrent);
  
  return elapsed_time;
}

bool AlmostEqual2sComplement(float A, float B, int maxUlps)
{
    // Make sure maxUlps is non-negative and small enough that the
    // default NAN won't compare as equal to anything.
    // assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
    int aInt = *(int*)&A;
    // Make aInt lexicographically ordered as a twos-complement int
    if (aInt < 0)
        aInt = 0x80000000 - aInt;
    // Make bInt lexicographically ordered as a twos-complement int
    int bInt = *(int*)&B;
    if (bInt < 0)
        bInt = 0x80000000 - bInt;
    int intDiff = abs(aInt - bInt);
    if (intDiff <= maxUlps)
        return true;
    return false;
}

cudaMemoryType getPtrLoc(void* ptr)
{
	cudaPointerAttributes ptrInfo;
	cudaMemoryType m;
	
	cudaPointerGetAttributes(&ptrInfo, ptr);
	
	m = ptrInfo.memoryType;
	if(ptrInfo.devicePointer!=ptr)
// 		Assume host
		m = cudaMemoryTypeHost;
	
	if((m!=cudaMemoryTypeHost) && (m!=cudaMemoryTypeDevice))
		//assume static memory is on host if no match
		m = cudaMemoryTypeHost;
	
	cudaGetLastError();
	return m;
}

void findBlockGrid(int npts, int blockSizeIn, size_t *block_size, size_t *grid_size)
{
	*block_size = blockSizeIn;
	*grid_size = npts/(*block_size);
	
// 	Temp fix for oversized grid
	while(*grid_size>65536)
	{
		*block_size *= 2;
		*grid_size = npts/(*block_size);
	}
	
	if(npts%*block_size)
		(*grid_size)++;
	
	return;
}

void findBlockGrid3(int npts, int maxThreadSize, dim3 *block_size, dim3 *grid_size)
{
	int numThreads;
	int numBlocks;
	int deviceID;
	cudaDeviceProp dp;

	cudaGetDevice(&deviceID);
	cudaGetDeviceProperties(&dp, deviceID);

	block_size->x = maxThreadSize;
	numBlocks = npts/block_size->x;
	
	// First increase number threads per block
	grid_size->x = numBlocks+1;
	while(grid_size->x>dp.maxGridSize[0] && (block_size->x<(maxThreadSize-32)))
	{
		block_size->x += 32;
		grid_size->x = ceil(npts/(1.0*block_size->x));
	}

	numBlocks = grid_size->x;
	grid_size->y = 1;
	while(grid_size->x>dp.maxGridSize[0])
	{
		grid_size->y++;
		grid_size->x = ceil(numBlocks/(1.0*grid_size->y));
	}
	
	/*
	if(npts%block_size->x)
		grid_size->x++;*/
	
	return;
}

void int3ToInts(int3 d, int *i, int N)
{
	i[0] = d.x;
	if(N>1)
	{
		i[1] = d.y;
		if(N>2)
			i[2] = d.z;
	}

	return;
}

int3 intsToInt3(int *i, int N, int def)
{
	int3 d;
	
	d.x = i[0];
	if(N>1)
		d.y = i[1] ;
	else
		d.y = def;
	
	if(N>2)
		d.z  = i[2];
	else
		d.z = def;

	return d;
}

int3 intsToInt3(int *i, int N)
{
	return intsToInt3(i,N,0);
}

float3 floatsToFloat3(float *i, int N)
{
	float3 d;
	
	d.x = i[0];
	if(N>1)
	{
		d.y = i[1] ;
		if(N>2)
			d.z  =i[2] ;
	}

	return d;
}

int bytes2mb(long m)
{
	m /= 1024;
	m /= 1024;
	
	return m;
}

long prodInt(int* x, int N)
{
	long p;
	int n;
	
	p = x[0];
	
	for(n=1; n<N; n++)
	 p *= x[n];
	
	return p;
}

size_t printAvailMem()
{
	return printAvailMem(0);
}

size_t printAvailMem(int quiet)
{
	size_t memFree, memTotal;
	
	cudaMemGetInfo(&memFree, &memTotal);
	if(!quiet)
		printf("%d/%d MB available on gpu\n", bytes2mb(memFree), bytes2mb(memTotal));
	
	return memFree;
}

int setGpu()
{
	int gpuID = 0;
	int numGpu;
	int n;
	size_t memAvail[4];
	size_t memMax = 0;	/* maximum available memory out of all gpus */

	cudaGetDeviceCount(&numGpu);

	for(n=0; n<numGpu; n++)
	{
		cudaSetDevice(n);
		memAvail[n] = printAvailMem(1); 
		if(memAvail[n]>memMax)
			gpuID = n;
	}

	return gpuID;
}

int findGpuMem(size_t mem)
{
	int gpuID = -1;
	int numGpu;
	int gpuIDOrig;	// Gpu function called from
	size_t memAvail;
	int n;
	
	cudaGetDevice(&gpuIDOrig);
	cudaGetDeviceCount(&numGpu);
	
	for(n=0; n<numGpu; n++)
	{
		cudaSetDevice(n);
		memAvail = printAvailMem(1); 
		if(memAvail>mem)
			gpuID = n;
		if(gpuID>=0)
			break;
	}
	
	cudaSetDevice(gpuIDOrig);
	
	
	return gpuID;
}

cudaMemcpyKind getCopyFlag(cudaMemoryType to, cudaMemoryType from)
{
	cudaMemcpyKind cf;
	
	if(from==cudaMemoryTypeDevice)
		if(to==cudaMemoryTypeHost)
			cf = cudaMemcpyDeviceToHost;
		else
			cf = cudaMemcpyDeviceToDevice;
	else if(from==cudaMemoryTypeHost)
		if(to==cudaMemoryTypeHost)
			cf = cudaMemcpyHostToHost;
		else
			cf = cudaMemcpyHostToDevice;
	else
	{
		fprintf(stderr, "Invalid data type %d\n", to);
		exit(EXIT_FAILURE);
	}
	
	return cf;
}

cudaMemcpyKind getCopyFlag(void* to, void* from)
{
	cudaMemoryType memFrom, memTo;
	memFrom = getPtrLoc(from);
	memTo = getPtrLoc(to);

	return getCopyFlag(memTo, memFrom);
}

int cuFree(void* ptr)
{
	int status = 0;

	cudaMemoryType mloc = getPtrLoc(ptr);

	if(mloc==cudaMemoryTypeHost)
#ifdef MATLAB_MEX_FILE
		mxFree(ptr);
#else
		cudaFreeHost(ptr);
#endif
	else if(mloc==cudaMemoryTypeDevice)
		cudaFree(ptr);
	else
	{
		fprintf(stderr, "Invalid memory location %d\n", mloc);
		status = 1;
	}

	if(!status)
		ptr = NULL;

	return status;
}

int cuMalloc(void** ptr, long size, cudaMemoryType mloc)
{
	int status = 0;

	if(mloc==cudaMemoryTypeHost)
#ifdef MATLAB_MEX_FILE
		*ptr = mxMalloc(size);
#else
		cudaMallocHost(ptr, size);
#endif
	else if(mloc==cudaMemoryTypeDevice)
		cudaMalloc(ptr, size);
	else
	{
		fprintf(stderr, "Invalid memory location %d\n", mloc);
		status = 1;
	}

	return status;
}

int cuCalloc(void** ptr, long size, cudaMemoryType mloc)
{
	int status;
	
	status = cuMalloc(ptr, size, mloc);

	if(!status)
		status = cuMemset(*ptr, 0, size);

	return status;
}

int cuMemset(void* ptr, int value, size_t size)
{
	cudaMemoryType mloc = getPtrLoc(ptr);
	int status = 0;

	if(mloc==cudaMemoryTypeHost)
		memset(ptr, value, size);
	else if(mloc==cudaMemoryTypeDevice)
		cudaMemset(ptr, value, size);
	else
	{
		fprintf(stderr, "Invalid memory location %d\n", mloc);
		status = 1;
	}

	return status;
}

__global__ void testCuComplex(cuComplex* x)
{
	unsigned int index = threadIdx.x + blockDim.x*blockIdx.x;

	return;
}

__global__ void testFloat(float * x)
{
	unsigned int index = threadIdx.x + blockDim.x*blockIdx.x;

	return;
}
