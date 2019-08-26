#include "convertendian.h"

#include <stdlib.h>

int inLittleEndian()
{
	/* The short value 1 has bytes (1,0) in little-endian mode and (0,1) in big-endian */
	short s = 1;
	return(((char*)&s)[0]) == 1;
}

int needEndianSwap(enum Endian endianDesired)
{
	int swap = 0;
	enum Endian endianSystem;

	if(inLittleEndian())
    endianSystem = LittleEndian;
	else
		endianSystem = BigEndian;

	swap = (endianDesired!=endianSystem);

	return swap;
}

void swapArrayEndian(void* buffer, int length, int ptsize)
{
	int i,j;
	unsigned char *tmp;
	unsigned char tmp2;
	int ptsizem1 = ptsize-1;

	tmp = (unsigned char *)buffer;

	for (i=0; i<length; i++)
	{
		for (j = 0; j < ptsize/2; j++)
		{
			tmp2 = tmp[ptsize*i +ptsizem1-j];
			tmp[ptsize*i +ptsizem1 - j] = tmp[ptsize*i +j];
			tmp[ptsize*i +j] = tmp2;
		 }
	}

	return;
}

void writeArray(void* array, size_t numPoints, int pointSize, FILE* file)
{
  if(array)
  {
    fwrite(&numPoints, sizeof(size_t), 1, file);
    fwrite(array, pointSize, numPoints, file);
  }
  else
  {
    numPoints = 0;
    fwrite(&numPoints, sizeof(size_t), 1, file);
  }
}

void readArray(void** array, int pointSize, FILE* file, enum Endian endian)
{
  size_t numPoints;
  fread(&numPoints, sizeof(size_t), 1, file);
  if(numPoints)
  {
    *array = malloc(numPoints*pointSize);
    fread(*array, pointSize, numPoints, file);
    if(needEndianSwap(endian))
      swapArrayEndian(*array, numPoints, pointSize);
  }
  else
    *array = NULL;
}
