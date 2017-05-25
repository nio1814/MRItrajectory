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
