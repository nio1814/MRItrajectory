#ifndef MP1UTIL
#define MP1UTIL

/**\defgroup	cudautil	CUDA*/

#include "cuComplex.h"

struct event_pair
{
  cudaEvent_t start;
  cudaEvent_t end;
};

/**
 \ingroup cudautil
 \brief	Check the GPU status and report last error
 \param[in]	kernel_name	Check point name
 \return	Last cuda error
*/
int check_launch(char * kernel_name);

/**
 \brief	Print message and exit
 \param[in]	msg	Exit message
*/
void printExit(char* msg);

void start_timer(struct event_pair * p);
float stop_timer(struct event_pair * p, char * kernel_name);

bool AlmostEqual2sComplement(float A, float B, int maxUlps);

/**
 \ingroup cudautil
 \brief	Determine if device or host pointer
 \param[in]	ptr	Pointer
 \return	Pointer location
*/
cudaMemoryType getPtrLoc(void* ptr);

/**
 \ingroup cudautil
 \brief	Determine the grid size 
 \param[in]	npts	Number of data points
 \param[in]	blockSize	Desired block size
 \param[out]	block_size	Allowable block size
 \param[out]	grid_size	Calculated grid size
*/
void findBlockGrid(int npts, int blockSize, size_t *block_size, size_t *grid_size);
void findBlockGrid3(int npts, int blockSizeIn, dim3 *block_size, dim3 *grid_size);

/**
 \ingroup	cudautil
 \brief	Convert int3 to an array of ints
 \param[in]	d	Original int3 variable
 \param[out]	i	Array of ints
 \param[in]	N	Number of used values
*/
/**
 \brief	Convert int3 to array of ints
*/
void int3ToInts(int3 d, int *i, int N);

/**
 \brief	Convert int array to int3
*/
int3 intsToInt3(int *i, int N);
int3 intsToInt3(int *i, int N, int def);

float3 floatsToFloat3(float *i, int N);

/**
 \brief Convert bytes to megabytes
 \param[in]	
*/
int bytes2mb(long m);

/**
 \brief	Calculate the product of an array of integers
*/
long prodInt(int* x, int N);

/**
 \brief	Display amount of available memory on GPU
 \return Memeory available (bytes)
*/
size_t printAvailMem(int quiet);
size_t printAvailMem();

/**
 \brief	Find GPU with adequate memory available
*/
int findGpuMem(size_t mem);

/**
 \param[in]	to	Memory locatin copying to
 \param[in]	to	Memory locatin copying from
*/
cudaMemcpyKind getCopyFlag(cudaMemoryType to, cudaMemoryType from);
cudaMemcpyKind getCopyFlag(void* to, void* from);

int cuMalloc(void** ptr, long size, cudaMemoryType mloc);
int cuCalloc(void** ptr, long size, cudaMemoryType mloc);
int cuFree(void* ptr);
int cuMemset(void* ptr, int value, size_t size);

__global__ void testCuComplex(cuComplex* x);
__global__ void testFloat(float * x);

#endif

