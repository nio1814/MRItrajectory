/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#ifndef ARRAYOPS_H
#define ARRAYOPS_H

void addfloat(float* arr, float addVal, int npts);
void addfloats(float *inputbuffer1, float *inputbuffer2, float *outputbuffer, long npts);
float sumfloats(float* arr, int npts);
void cumulativeSum(float* arr, float* arrOut, int numPts);

void scalefloats(float *data, long npts, float multiplier);
void multiplyfloats(float* arr1, float* arr2, float* arrOut, int numPts);

void subtractArray(float* arr1, float* arr2, float* arrOut, int numPts);

void diffArray(float* arr, float* arrOut, int numPts);

void interpolatefloats(float *arrx, float *arry, int nptsIn, float *x, float *y, int nptsOut);

/**
 \brief	Set all values to a given value
  @param[in]	data		data values
  @param[in]	npts		Number of points
  \param[in]	value	Value to set to
 */
void setValues(float *data, long npts, float value);

/** L2 norm
\ingroup	Math
\param[in]	arr	Input vector
\param[in]	npts	Length of input vector
\return	L2 norm
*/
float norm2(float* arr, int npts);

/**
 \brief	Find the minimum value in a float array
 \param[in]	arr1	Array of floats
 \param[in]	numPts	Array length
 \return	minimum value in array
*/
float minValue(const float* arr1, int numPts);

/**
\brief	Calculate the maximum value of the data array.
  @param[in]	data		data values i
  @param[in]	npts		Number of points.
  \return maximum of the array.
*/
float maxValue(const float *data, long npts);

/**
\ingroup arrayops
  \brief	Calculate the maximum absolute value of the data array.
  \param[in]	rdata	real data values i
  \param[in]	idata	imag data values i
  \param[in]	npts	Number of points.
  \return maximum of the array.
*/
float maxMagnitude(float *rdata, float *idata, long npts);

float distance(float* x, float* y, int npts);

float dot(float* x, float *y, int npts);

void swapValues(float *array, int index1, int index2);

void scalecomplex(float *rdata, float *idata, float rmult, float imult, long npts);

#endif // ARRAYOPS_H
