/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#include "gridding.h"

extern "C" {
#include "trajectory.h"
#include "numericalrecipes.h"
#include "arrayops.h"
}

Gridding::Gridding(const Trajectory *trajectory) :
	m_trajectory(trajectory), m_oversamplingFactor(1.5),
	m_kernelWidth(4)
{
	float beta = M_PI*sqrt(powf(m_kernelWidth/m_oversamplingFactor,2)*powf(m_oversamplingFactor-0.5,2)-0.8);

	float kb1 = -0.4491*m_oversamplingFactor - .0298;
	float kb2 = 0.2522*m_oversamplingFactor - .0378;
	float maxError = pow(10.0, kb1*m_kernelWidth+kb2);

	int numKernelPoints = sqrt(0.37/maxError)/m_oversamplingFactor;
	int numKernelPointsHalf = ceil(0.5*numKernelPoints*m_kernelWidth);

	float scale = 1/besseli(0, beta);
	float kernelSum = 0;
	for(int n=0; n<numKernelPointsHalf; n++)
	{
		float x = n/(float)numKernelPoints;
		float value = besseli(0, beta*sqrt(1-pow(2*x/m_kernelWidth,2)))*scale;
		kernelSum += 2*value;
		m_kernelLookupTable.push_back(value);
	}
	kernelSum -= m_kernelLookupTable[0];
	kernelSum /= numKernelPoints*sqrt(m_oversamplingFactor);
	scalefloats(m_kernelLookupTable.data(), m_kernelLookupTable.size(), 1/kernelSum);


	m_deapodization.clear();
	scale = -1/(2*beta)*(exp(-beta)-exp(beta));

	float minSpatialResolution = INFINITY;
	for(int d=0; d<numDimensions(); d++)
		minSpatialResolution = std::min(minSpatialResolution, m_trajectory->spatialResolution[d]);

	for(int d=0; d<numDimensions(); d++)
	{
		int dimension = m_oversamplingFactor*m_trajectory->imageDimensions[d];
		dimension += dimension%2;
		m_gridDimensions.push_back(dimension);
		m_normalizedToGridScale.push_back(m_trajectory->spatialResolution[d]/minSpatialResolution*m_gridDimensions[d]);
	}

	m_deapodization.resize(numDimensions());
	for(int d=0; d<m_trajectory->dimensions; d++)
	{
		for(int n=0; n<m_gridDimensions[d]; n++)
		{
			float arg = M_PI*m_kernelWidth*(-m_gridDimensions[d]/2.0f + n)/m_gridDimensions[d];
			arg *= arg;
			arg -= beta*beta;

			float value;
			if(arg<0)
				value = -1.0f/(2.0f*sqrtf(-arg))*(expf(-sqrtf(-arg))-expf(sqrtf(-arg)));
			else
				value = sinf(sqrtf(arg))/sqrtf(arg);
			value = scale/value;
			m_deapodization[d].push_back(value);
		}
	}
}

std::vector<int> Gridding::imageDimensions()
{
	std::vector<int> imageDims;
	for(int d=0; d<m_trajectory->dimensions; d++)
		imageDims.push_back(m_trajectory->imageDimensions[d]);

	return imageDims;
}

void Gridding::nearestGriddedPoint(const std::vector<float> &ungriddedPoint, std::vector<float> &griddedPoint)
{
	for(int d=0; d<numDimensions(); d++)
		griddedPoint[d] = roundf(m_normalizedToGridScale[d]*ungriddedPoint[d])/m_normalizedToGridScale[d];
}

float Gridding::lookupKernelValue(float x)
{
	float kernelIndexDecimal = fabsf((x)/(float)m_kernelWidth/2*m_kernelLookupTable.size());
	int kernelIndex = (int)kernelIndexDecimal;
	float secondPointFraction = kernelIndexDecimal - kernelIndex;

	return m_kernelLookupTable[kernelIndex]*(1-secondPointFraction) + m_kernelLookupTable[kernelIndex+1]*secondPointFraction;
}

std::vector<std::vector<float> > Gridding::kernelNeighborhood(const std::vector<float> &ungriddedPoint, const std::vector<float> &griddedPoint)
{
	int lookupPoints = m_kernelLookupTable.size();
	std::vector<std::vector<float> > kernelValues(numDimensions());

	for(int d=0; d<numDimensions(); d++)
	{
		float offset = (griddedPoint[d]-ungriddedPoint[d])*m_normalizedToGridScale[d];
		for(int n=0; n<m_kernelWidth; n++)
		{
			int m = n-m_kernelWidth/2+1;
			float kernelIndexDecimal = fabsf((m+offset)/(float)m_kernelWidth/2*lookupPoints);
			int kernelIndex = (int)kernelIndexDecimal;
			float secondPointFraction = kernelIndexDecimal - kernelIndex;
			kernelValues[d].push_back(m_kernelLookupTable[kernelIndex]*(1-secondPointFraction) + m_kernelLookupTable[kernelIndex+1]*secondPointFraction);
		}
	}

	return kernelValues;
}

size_t multiToSingleIndex(const std::vector<int> &indices, const std::vector<int>& size)
{
	long index = indices[0];
	long lowerDimensionSize = size[0];
	for(size_t n=1; n<indices.size(); n++)
	{
		index += indices[n]*lowerDimensionSize;
		lowerDimensionSize *= size[n];
	}

	return index;
}

void singleToMultiIndex(long index, const std::vector<int>& size, std::vector<int>& indices)
{
	int lowerDimensions = 1;
	for(size_t n=0; n<size.size(); n++)
	{
		indices[n] = (index/lowerDimensions)%size[n];
		lowerDimensions *= size[n];
	}
}

MRdata *Gridding::grid(const MRdata &ungriddedData)
{
	std::vector<float> ungriddedPoint(3);
	std::vector<float> griddedPoint(3);

	int dimensionStart[3] = {0,0,0};
	int dimensionEnd[3] = {1,1,1};
	int offset[3] = {0,0,0};

	MRdata* griddedData = new MRdata(m_gridDimensions, numDimensions());

	std::vector<int> gridDimensions3 = m_gridDimensions;
	if(numDimensions()==2)
		gridDimensions3.push_back(1);
	std::vector<float> one;
	one.push_back(1);

	for(size_t n=0; n<ungriddedData.points(); n++)
	{
		complexFloat ungriddedDataValue = ungriddedData.signalValue(n);
		int readoutPoint = n%m_trajectory->readoutPoints;
		int readout = n/m_trajectory->readoutPoints;

		trajectoryCoordinates(readoutPoint, readout, m_trajectory, ungriddedPoint.data());

		nearestGriddedPoint(ungriddedPoint, griddedPoint);

		for(int d=0; d<numDimensions(); d++)
		{
			int gridPointCenter = (int)((griddedPoint[d]+.5)*m_normalizedToGridScale[d]);
			int start = gridPointCenter-m_kernelWidth/2+1;
			dimensionStart[d] = std::max(0, start);
			offset[d] = dimensionStart[d]-start;
			dimensionEnd[d] = std::min(gridPointCenter+m_kernelWidth/2, m_gridDimensions[d]);
		}

		std::vector<std::vector<float> > kernelValues = kernelNeighborhood(ungriddedPoint, griddedPoint);
		std::vector<int> gridIndex(3);
		for(int d=numDimensions(); d<3; d++)
			kernelValues.push_back(one);

		for(int gz=dimensionStart[2]; gz<dimensionEnd[2]; gz++)
		{
			gridIndex[2] = gz;
			float kernelZ = kernelValues[2][gz-dimensionStart[2]+offset[2]];
			for(int gy=dimensionStart[1]; gy<dimensionEnd[1]; gy++)
			{
				gridIndex[1] = gy;
				float kernelY = kernelValues[1][gy-dimensionStart[1]+offset[1]];
				for(int gx=dimensionStart[0]; gx<dimensionEnd[0]; gx++)
				{
					gridIndex[0] = gx;
					float kernelX = kernelValues[0][gx-dimensionStart[0]+offset[0]];
					long griddedDataIndex = multiToSingleIndex(gridIndex, gridDimensions3);
					griddedData->setSignalValue(griddedDataIndex, kernelX*kernelY*kernelZ*ungriddedDataValue);
				}
			}
		}
	}

	return griddedData;
}

void Gridding::deapodize(MRdata &oversampledImage)
{
	complexFloat* signal = oversampledImage.signalPointer();

	std::vector<int> gridIndex(numDimensions());
	for(size_t n=0; n<oversampledImage.points(); n++)
	{
		singleToMultiIndex(n, m_gridDimensions, gridIndex);
		for(int d=0; d<numDimensions(); d++)
			signal[n] *= m_deapodization[d][gridIndex[d]];
	}
}

MRdata* Gridding::kSpaceToImage(const MRdata &ungriddedData)
{
	MRdata* image = grid(ungriddedData);
	image->fftShift();
	image->fft(FFTW_BACKWARD);
	image->fftShift();

	deapodize(*image);
	image->crop(imageDimensions());
	image->writeToOctave("temp");

	return image;
}

MRdata *Gridding::imageToKspace(const MRdata &image)
{
	MRdata griddedKspace = image;
	griddedKspace.pad(m_gridDimensions);
	deapodize(griddedKspace);

	griddedKspace.fftShift();
	griddedKspace.fft(FFTW_FORWARD);
	griddedKspace.fftShift();

	MRdata* ungriddedKspace = inverseGrid(griddedKspace);

	return ungriddedKspace;
}

int Gridding::numDimensions()
{
	return m_trajectory->dimensions;
}

std::vector<int> Gridding::gridDimensions()
{
	return m_gridDimensions;
}


