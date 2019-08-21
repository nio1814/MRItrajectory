/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#include "arrayops.h"

#include <math.h>
#include <stdio.h>

void addfloat(float* arr, float addVal, int npts)
{
	int i;

	for(i=0; i<npts; i++)
		arr[i] += addVal;

	return;
}

void addfloats(float *inputbuffer1, float *inputbuffer2, float *outputbuffer, long npts)
{
	long count;
	for (count=0; count < npts; count++)
		*outputbuffer++ = (*inputbuffer1++) + (*inputbuffer2++);
}

float sumfloats(float* arr, int npts)
{
	float sumVal;
	int i;

	sumVal = 0;

	for(i=0; i<npts; i++)
		sumVal += arr[i];

	return sumVal;
}

void cumulativeSum(float* arr, float* arrOut, int numPts)
{
	int i;
	float sum = 0.0f;

	for(i=0; i<numPts; i++)
	{
		sum += arr[i];
		arrOut[i] = sum;
	}

	return;
}

void scalefloats(float *data, long npts, float multiplier)
{
	long count;

	for (count=0; count < npts; count++)
	{
		(*data) *= multiplier;
		data++;
	}

	return;
}

void multiplyfloats(float* arr1, float* arr2, float* arrOut, int numPts)
{
	int i;

	for(i=0; i<numPts; i++)
		arrOut[i] = arr1[i]*arr2[i];

	return;
}

void subtractArray(float* arr1, float* arr2, float* arrOut, int numPts)
{
	int i;

	for(i=0; i<numPts; i++)
		arrOut[i] = arr1[i]-arr2[i];

	return;
}

void diffArray(float* arr, float* arrOut, int numPts)
{
	int i;

	for(i=0; i<(numPts-1); i++)
		arrOut[i] = arr[i+1]-arr[i];

	return;
}

float norm2(float* arr, int npts)
{
	float result=0;
	int i;
	for(i=0; i<npts; i++)
		result += arr[i]*arr[i];
	result = sqrt(result);
	return result;
}

float distance(float* x, float* y, int npts)
{
	float result=0;
	int i;

	for(i=0; i<npts; i++)
		result += (y[i]-x[i])*(y[i]-x[i]);
	result = sqrt(result);

	return result;
}

float dot(float* x, float *y, int npts)
{
	int i;
	float result = 0;

	for(i=0; i<npts; i++)
		result += x[i]*y[i];

	return result;
}

void swapValues(float *array, int index1, int index2)
{
	float value1 = array[index1];
	array[index1] = array[index2];
	array[index2] = value1;
}

void scalecomplex(float *rdata, float *idata, float rmult, float imult, long npts)
{
	long count;
	float rval, ival;

	for (count=0; count < npts; count++)
	{
		rval = (*rdata)*rmult - (*idata)*imult;
		ival = (*rdata)*imult + (*idata)*rmult;
		*rdata++ = rval;
		*idata++ = ival;
	}

	return;
}

float interpfloat(float *arrx, float *arry, int npts, float x)
{
	int n;
	int nLo=0, nHi=npts-1;

	/* Check for monotonicity */
	if(arrx[1]>arrx[0])
	{
		for(n=1; n<npts-1; n++)
			if(arrx[n+1]<=arrx[n])
			{
				fprintf(stderr, "Error: interpfloats() x is not monotonically increasesing");
				return 0;
			}
	}
	else if(arrx[1]<arrx[0])
  {
		for(n=1; n<npts-1; n++)
			if(arrx[n+1]>=arrx[n])
			{
				fprintf(stderr, "Error: interpfloats() x is not monotonically decreasing");
				return 0;
			}
  }
  else
	{
		fprintf(stderr, "Error: interpfloats() x is not monotonically decreasing or increasesing");
		return 0;
	}

	do
	{
		n = .5*(nLo+nHi);
		if(x>arrx[n])
			nLo = n;
		else
			nHi = n;
	} while((nHi-nLo)>1);


	return ((arrx[nHi]-x)*arry[nLo]+(x-arrx[nLo])*arry[nHi])/(arrx[nHi]-arrx[nLo]);
}


void interpolatefloats(float *arrx, float *arry, int nptsIn, float *x, float *y, int nptsOut)
{
	int n;
	for(n=0; n<nptsOut; n++)
		y[n] = interpfloat(arrx, arry, nptsIn, x[n]);
}

float minValue(const float *arr1, int numPts)
{
  float mv;
  int i;

  mv = arr1[0];

  for(i=1; i<numPts; i++)
  {
    if(arr1[i]<mv)
      mv = arr1[i];
  }

  return mv;
}

float maxValue(const float *data, long npts)
{
  float maxval;
  long count;

  maxval = *data;
  for (count=0; count < npts; count++)
  {
    if (*data > maxval)
      maxval = *data;
    data++;
  }

  return (maxval);
}

float maxMagnitude(float *rdata, float *idata, long npts)
{
  float valsq, maxval,maxvalsq;
  long count;

  maxvalsq = 0;
  for (count=0; count < npts; count++)
    {
    valsq = *rdata * *rdata + *idata * *idata;
    if (valsq > maxvalsq)
      maxvalsq = valsq;
    rdata++;
    idata++;
    }
  maxval = sqrt(maxvalsq);

  return (maxval);
}

void setValues(float *data, long npts, float value)
{
  long count;
  for (count=0; count < npts; count++)
    *data++ = value;
}
