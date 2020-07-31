/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#include "mrdata.h"

extern "C"
{
#include "arrayops.h"
}

#include <string.h>

size_t dimensionsToPoints(std::vector<int> dimensions)
{
	size_t points = 1;
	for(size_t n=0; n<dimensions.size(); n++)
		points *= dimensions[n];

	return points;
}

MRdata::MRdata(const MRdata *from)
{
	*this = *from;
}

MRdata::MRdata(const std::vector<int> dimensions, const int numImagingDimensions, const std::vector<complexFloat> &data) :
  m_dimensions(dimensions), m_numImagingDimensions(numImagingDimensions)
{
	if(data.empty())
		m_signal.resize(points());
	else
    m_signal = data;
}

//MRdata::MRdata(PyArrayObject *data)
//{
//  initializeNumpy();

////  PyArrayObject* array = reinterpret_cast<PyArrayObject*>(PyArray_Cast(reinterpret_cast<PyArrayObject*>(data), NPY_COMPLEX64));
//  PyArrayObject* array = reinterpret_cast<PyArrayObject*>(PyArray_Cast(data, NPY_COMPLEX64));
//  const npy_intp *dimensions = PyArray_DIMS(array);
//  int numPoints = 1;

//  for (int d=0; d<PyArray_NDIM(array); d++) {
//    m_dimensions.push_back(dimensions[d]);
//    numPoints *= dimensions[d];
//  }

//  m_signal.reserve(numPoints);
//  complexFloat* arrayData = reinterpret_cast<complexFloat*>(PyArray_DATA(array));
//  std::copy(arrayData, arrayData + numPoints, std::back_inserter(m_signal));

//  m_numImagingDimensions = m_dimensions.size() > 3 ? 3 : m_dimensions.size();
//}

size_t MRdata::points() const
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

void MRdata::addSignalValue(size_t index, const complexFloat &value)
{
	m_signal[index] += value;
}

complexFloat MRdata::signalValue(int index) const
{
	return m_signal[index];
}

complexFloat *MRdata::signalPointer()
{
	return m_signal.data();
}

std::vector<complexFloat> MRdata::signal()
{
	return m_signal;
}

std::vector<int> MRdata::dimensions()
{
	return m_dimensions;
}

void MRdata::fftShift()
{
	int n;

	int index;
	int indexOriginal[3] = {0,0,0};
	int shift[3] = {0,0,0};

  int numShifts = static_cast<int>(points()) / m_dimensions[0];

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
			m_signal[indexShift+m] = m_signal[index+m];
			m_signal[index+m] = temp;
		}
	}
}

FFTplan MRdata::planFFT(int direction, std::vector<int> outputDimensions)
{
	bool inPlace = true;
	for(int n=0; n<m_numImagingDimensions; n++)
		if(m_dimensions[n]!=outputDimensions[n])
			inPlace = false;

  const size_t pointsOutput = dimensionsToPoints(outputDimensions);

	fftwf_iodim transformStrides[1];
	transformStrides[0].n = 1;
  transformStrides[0].is = static_cast<int>(points());
  transformStrides[0].os = static_cast<int>(pointsOutput);

	fftwf_iodim transformDimensions[3];
	int n = m_numImagingDimensions-1;
	transformDimensions[n].n = m_dimensions[m_numImagingDimensions-n-1];
	transformDimensions[n].is = 1;
	transformDimensions[n].os = 1;

	for(int n=m_numImagingDimensions-2; n>=0; n--)
	{
		transformDimensions[n].n = m_dimensions[m_numImagingDimensions-n-1];
		transformDimensions[n].is = m_dimensions[n]*transformDimensions[n+1].is;
		transformDimensions[n].os = outputDimensions[n]*transformDimensions[n+1].os;
	}

//	std::vector<complexFloat> signalIn(points());
//	std::vector<complexFloat> signalOut;
	fftwf_complex* signalIn = (fftwf_complex*)fftwf_malloc(points()*sizeof(fftwf_complex));
	fftwf_complex* signalOut;
	if(inPlace)
//		signalOutPointer = signalIn.data();
		signalOut = signalIn;
	else
	{
//		signalOut.resize(pointsOutput);
//		signalOutPointer = signalOut.data();
		signalOut = (fftwf_complex*)fftwf_malloc(pointsOutput*sizeof(fftwf_complex));
	}

	FFTplan plan;
//	plan.plan = fftwf_plan_guru_dft(m_numImagingDimensions, transformDimensions, 1, transformStrides, (fftwf_complex*)signalIn.data(), (fftwf_complex*)signalOutPointer, direction, FFTW_ESTIMATE);
	plan.plan = fftwf_plan_guru_dft(m_numImagingDimensions, transformDimensions, 1, transformStrides, signalIn, signalOut, direction, FFTW_ESTIMATE);
	plan.inPlace = inPlace;

	fftwf_free(signalIn);
	if(!inPlace)
		fftwf_free(signalOut);

	return plan;
}

void MRdata::fft(int direction, std::vector<int> outputDimensions)
{
	if(outputDimensions.empty())
		outputDimensions = m_dimensions;
	FFTplan plan = planFFT(direction, outputDimensions);

	size_t pointsOutput = dimensionsToPoints(outputDimensions);

	/*std::vector<complexFloat> signalNew;
	complexFloat* signalNewPointer = NULL;
	if(plan.inPlace)
		signalNewPointer = m_signal.data();
	else
	{
		signalNew.resize(pointsOutput);
		signalNewPointer = signalNew.data();
	}

	fftwf_execute_dft(plan.plan, (fftwf_complex*)m_signal.data(), (fftwf_complex*)signalNewPointer);

	if(!plan.inPlace)
		m_signal = signalNew;
*/

	fftwf_complex* signalIn = (fftwf_complex*)fftwf_malloc(points()*sizeof(fftwf_complex));
	fftwf_complex* signalOut;
	if(plan.inPlace)
		signalOut = signalIn;
	else
	{
		signalOut = (fftwf_complex*)fftwf_malloc(pointsOutput*sizeof(fftwf_complex));
	}

	memcpy(signalIn, m_signal.data(), points()*sizeof(fftwf_complex));
	fftwf_execute_dft(plan.plan, signalIn, signalOut);

  memcpy(m_signal.data(), reinterpret_cast<complexFloat*>(signalOut), pointsOutput*sizeof(fftwf_complex));

	fftwf_free(signalIn);
	if(!plan.inPlace)
		fftwf_free(signalOut);

	m_dimensions = outputDimensions;

  scalefloats(reinterpret_cast<float*>(m_signal.data()), 2*points(), 1.0 / std::sqrt(points()));
}

void MRdata::crop(std::vector<int> newSize)
{
	int pointsCrop = 1;
	for(size_t n=0; n<newSize.size(); n++)
		pointsCrop *= newSize[n];

	std::vector<complexFloat> signalCropped(pointsCrop);

	int numCopies = pointsCrop/newSize[0];

	int multiIndexCrop[3] = {0,0,0};
	int indexOriginal[3] = {0,0,0};
	for(int n=0; n<numCopies; n++)
	{
		multiIndexCrop[1] = n%newSize[1];
		multiIndexCrop[2] = n/newSize[1];

		for(int s=0; s<m_numImagingDimensions; s++)
    {
			indexOriginal[s] = multiIndexCrop[s] - (newSize[s]-m_dimensions[s])/2;

			int index = (indexOriginal[2]*m_dimensions[1] + indexOriginal[1])*m_dimensions[0] + indexOriginal[0];
			int indexCrop = (multiIndexCrop[2]*newSize[1] + multiIndexCrop[1])*newSize[0] + multiIndexCrop[0];

			for(int m=0; m<newSize[0]; m++)
				signalCropped[indexCrop+m] = m_signal[index+m];
    }
  }

	m_signal = signalCropped;
	for(int d=0; d<m_numImagingDimensions; d++)
		m_dimensions[d] = newSize[d];
}

void MRdata::pad(std::vector<int> newSize)
{
	/* Calculate new total number of points */
	size_t pointsPad = 1;
	int shift[3] = {0,0,0};
	bool doPad = false;	// perform padding
	for(int d=0; d<m_numImagingDimensions; d++)
	{
    doPad |= m_dimensions[d] < newSize[d];
		pointsPad *= newSize[d];
    shift[d] = static_cast<int>(std::ceil((newSize[d] - m_dimensions[d])/2.0));
	}

	if(!doPad)
		return;

	std::vector<int> dimensionsOriginal = m_dimensions;
	for(int d=m_numImagingDimensions; d<3; d++)
	{
		newSize.push_back(1);
		dimensionsOriginal.push_back(1);
	}

	/* Make new signal */
	std::vector<complexFloat> signalPad = std::vector<complexFloat>(pointsPad);

	int nx = shift[0];

	for(int nyOld=0; nyOld<dimensionsOriginal[1]; nyOld++)
	{
		int ny = nyOld + shift[1];
			for(int nzOld=0; nzOld<dimensionsOriginal[2]; nzOld++)
			{
				int nz = nzOld + shift[2];

				size_t n = (nz*newSize[1] + ny)*newSize[0] + nx;
				int nOld = (nzOld*dimensionsOriginal[1] + nyOld)*dimensionsOriginal[0];

				for(int m=0; m<dimensionsOriginal[0]; m++)
//				blas_ccopy(b_N, &signal[nOld+s*npts], b_inc, &signalPad[n+s*pointsPad], b_inc);
					signalPad[n+m] = m_signal[nOld+m];
			}
	}

	m_signal = signalPad;
	m_dimensions = newSize;
}

bool MRdata::writeToOctave(std::string filename) const
{
	FILE* file = fopen(filename.c_str(), "w");
	if(!file)
	{
		fprintf(stderr, "Failed to open %s", filename.c_str());
		return true;
	}

	fprintf(file, "# Phantom Test\n");
	int periodPosition = filename.find('.');
	std::string name;
	if(periodPosition>0)
		name = filename.substr(0,periodPosition);
	else
		name = filename;
	fprintf(file, "# name: %s\n", name.c_str());
	fprintf(file, "# type: complex matrix\n");
	fprintf(file, "# ndims: %zu\n", m_dimensions.size());
	for(size_t d=0; d<m_dimensions.size(); d++)
		fprintf(file, " %d", m_dimensions[d]);
	fprintf(file, "\n");
//	QString line = QString(" %1 %2 %3\n").arg(imageSize[0]).arg(imageSize[1]).arg(imageSize[2]);

	for (size_t n=0; n<m_signal.size(); n++)
		fprintf(file, " (%f,%f)\n", m_signal[n].real(), m_signal[n].imag());

	fclose(file);

	return false;
}

