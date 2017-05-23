/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#include "mrdata.h"

MRdata::MRdata(std::vector<int> dimensions, int numImagingDimensions, const std::vector<complexFloat> &data) : m_numImagingDimensions(numImagingDimensions)
{
	if(data.empty())
		m_signal.resize(points());
	else
		m_signal = data;
}

long MRdata::points() const
{
	long N = 0;
	for(int d=0; d<m_numImagingDimensions; d++)
	{
		if(!N)
			N = 1;
		N *= m_dimensions[d];
	}

	return N;
}

void MRdata::setSignalValue(int index, const complexFloat &value)
{
	m_signal[index]	= value;
}

complexFloat MRdata::signalValue(int index) const
{
	return m_signal[index];
}

complexFloat *MRdata::signal()
{
	return m_signal.data();
}

void MRdata::fftShift()
{
	int n;

	int index;
	int indexOriginal[3] = {0,0,0};
	int shift[3] = {0,0,0};

	int numShifts = points()/m_dimensions[0];

	for(n=0; n<numShifts; n++)
	{
		indexOriginal[1] = n%m_dimensions[1];
		indexOriginal[2] = n/m_dimensions[1];

		for(int d=0; d<m_numImagingDimensions; d++)
			shift[d] = (indexOriginal[d]+m_dimensions[d]/2)%m_dimensions[d];

		index = (indexOriginal[2]*m_dimensions[1]+indexOriginal[1])*m_dimensions[0] + indexOriginal[0];
		int indexShift = (shift[2]*m_dimensions[1]+shift[1])*m_dimensions[0] + shift[0];

		for(int m=0; m<m_dimensions[0]/2; m++)
		{
			complexFloat temp = m_signal[indexShift+m];
			m_signal[indexShift+m] = m_signal[indexShift];
			m_signal[index] = temp;
		}
//		blas_ccopy(b_N, &signal[indexShift+setshift], 1, temp, 1);
//		blas_ccopy(b_N, &signal[index+setshift], 1, &signal[indexShift+setshift], 1);
//		blas_ccopy(b_N, temp, 1, &signal[index+setshift], 1);
	}
}

void MRdata::crop(std::vector<int> newSize)
{
	int n;

	int numCopies;
	int multiIndexCrop[3] = {0,0,0};
	int indexOriginal[3] = {0,0,0};
	int s;
	int index;

	int b_N = newSize[0];

	int pointsCrop = 1;
	for(n=0; n<m_numImagingDimensions; n++)
		pointsCrop *= newSize[n];

	std::vector<complexFloat> signalCropped;

	numCopies = pointsCrop/newSize[0];

	for(n=0; n<numCopies; n++)
	{
		multiIndexCrop[1] = n%newSize[1];
		multiIndexCrop[2] = n/newSize[1];

		for(s=0; s<m_numImagingDimensions; s++)
			indexOriginal[s] = multiIndexCrop[s] - (newSize[s]-m_dimensions[s])/2;

			index = (indexOriginal[2]*m_dimensions[1] + indexOriginal[1])*m_dimensions[0] + indexOriginal[0];
			int indexCrop = (multiIndexCrop[2]*newSize[1] + multiIndexCrop[1])*newSize[0] + multiIndexCrop[0];

			for(int m=0; m<newSize[0]; m++)
				signalCropped[indexCrop] = m_signal[index];
	}

	m_signal = signalCropped;
	for(int d=0; d<m_numImagingDimensions; d++)
		m_dimensions[d] = newSize[d];
}

