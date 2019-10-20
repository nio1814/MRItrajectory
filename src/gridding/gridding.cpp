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

#include <cmath>
#include <algorithm>


Gridding::Gridding(const Trajectory *trajectory, float oversamplingRatio, int kernelWidth) :
	m_trajectory(trajectory), m_oversamplingFactor(oversamplingRatio),
	m_kernelWidth(kernelWidth)
{
	float beta = M_PI*sqrt(powf(m_kernelWidth/m_oversamplingFactor,2)*powf(m_oversamplingFactor-0.5,2)-0.8);


	int numKernelPoints = 30*m_kernelWidth;

	float scale = 1/besseli(0, beta);
	float kernelSum = 0;
	for(int n=0; n<=numKernelPoints; n++)
	{
		float x = 2*n/(float)numKernelPoints;
		float value = besseli(0, beta*sqrt(1-pow(2*x/m_kernelWidth,2)))*scale;
		kernelSum += 2*value;
		m_kernelLookupTable.push_back(value);
	}
	kernelSum -= m_kernelLookupTable[0];
	kernelSum /= .5*numKernelPoints*sqrt(m_oversamplingFactor);
	scalefloats(m_kernelLookupTable.data(), m_kernelLookupTable.size()-1, 1/kernelSum);

	m_deapodization.clear();
	scale = -1/(2*beta)*(exp(-beta)-exp(beta));

	float minSpatialResolution = INFINITY;
	for(int d=0; d<numDimensions(); d++)
		minSpatialResolution = std::min(minSpatialResolution, m_trajectory->spatialResolution[d]);
	m_coordinateScale = .5*minSpatialResolution/5;
	for(int d=0; d<numDimensions(); d++)
	{
    int dimension = static_cast<int>(std::roundf(m_oversamplingFactor*m_trajectory->imageDimensions[d]));
		dimension += dimension%2;
		m_gridDimensions.push_back(dimension);
		m_normalizedToGridScale.push_back(m_trajectory->spatialResolution[d]/minSpatialResolution*m_gridDimensions[d]);
	}

	m_deapodization.resize(numDimensions());
  for(int d=0; d<m_trajectory->numDimensions; d++)
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
  for(int d=0; d<m_trajectory->numDimensions; d++)
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
	int lookupPoints = m_kernelLookupTable.size()-1;
	std::vector<std::vector<float> > kernelValues(numDimensions());

	for(int d=0; d<numDimensions(); d++)
	{
		float offset = (griddedPoint[d]-ungriddedPoint[d])*m_normalizedToGridScale[d]-1;
		offset += ungriddedPoint[d]>=griddedPoint[d] ? 1 : 0;
		for(int n=0; n<m_kernelWidth; n++)
		{
			int m = n-m_kernelWidth/2+1;
			float kernelOffset = m+offset;
			float kernelIndexDecimal = fabsf(kernelOffset/(float)m_kernelWidth*2*(lookupPoints-1));
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

MRdata *Gridding::grid(MRdata &inputData, Direction direction)
{
	std::vector<float> ungriddedPoint(3);
	std::vector<float> griddedPoint(3);

	int dimensionStart[3] = {0,0,0};
	int dimensionEnd[3] = {1,1,1};
	int offset[3] = {0,0,0};

	MRdata* ungriddedData;
	MRdata* griddedData;

	if(direction==Forward)
	{
		ungriddedData = &inputData;
		griddedData = new MRdata(m_gridDimensions, numDimensions());
	}
	else
	{
		std::vector<int> trajectorySize(2);
    trajectorySize[0] = m_trajectory->numReadoutPoints;
    trajectorySize[1] = m_trajectory->numReadouts;
		ungriddedData = new MRdata(trajectorySize, numDimensions());
		griddedData = &inputData;
	}

	std::vector<int> gridDimensions3 = m_gridDimensions;
	if(numDimensions()==2)
		gridDimensions3.push_back(1);
	std::vector<float> one;
	one.push_back(1);

	for(size_t n=0; n<ungriddedData->points(); n++)
	{
		complexFloat ungriddedDataValue;
		if(direction==Forward)
			ungriddedDataValue = ungriddedData->signalValue(n);
		else
			ungriddedDataValue = 0;
    int readoutPoint = n%m_trajectory->numReadoutPoints;
    int readout = n/m_trajectory->numReadoutPoints;

		float densityCompensation;
		trajectoryCoordinates(readoutPoint, readout, m_trajectory, ungriddedPoint.data(), &densityCompensation);
		scalefloats(ungriddedPoint.data(), numDimensions(), m_coordinateScale);
		ungriddedDataValue *= densityCompensation;

		nearestGriddedPoint(ungriddedPoint, griddedPoint);

		for(int d=0; d<numDimensions(); d++)
		{
      const int gridPointCenter = static_cast<int>(std::round((griddedPoint[d]+.5)*m_normalizedToGridScale[d]));
			int start = gridPointCenter-m_kernelWidth/2;
			start += ungriddedPoint[d]>=griddedPoint[d] ? 1 : 0;
			dimensionStart[d] = std::max(0, start);
			offset[d] = dimensionStart[d]-start;
			dimensionEnd[d] = std::min(start+m_kernelWidth, m_gridDimensions[d]);
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
					if(direction==Forward)
						griddedData->addSignalValue(griddedDataIndex, kernelX*kernelY*kernelZ*ungriddedDataValue);
					else
						ungriddedDataValue += kernelX*kernelY*kernelZ*griddedData->signalValue(griddedDataIndex);
				}
			}
		}
		if(direction==Inverse)
			ungriddedData->setSignalValue(n, ungriddedDataValue);
	}

	if(direction==Forward)
		return griddedData;
	else
		return ungriddedData;
}

MRdata *Gridding::conjugatePhaseForward(const MRdata &ungriddedData)
{
  std::vector<int> imageDimensions(m_trajectory->imageDimensions, &m_trajectory->imageDimensions[m_trajectory->numDimensions]);
	MRdata* griddedData = new MRdata(imageDimensions, numDimensions());
	float imageScale = 1.0f/std::sqrt(griddedData->points());

	std::vector<int> imageCenter;
	std::vector<float> axisScale;
	for(int d=0; d<numDimensions(); d++)
	{
		axisScale.push_back(m_normalizedToGridScale[d]/m_gridDimensions[d]*m_coordinateScale);
		imageCenter.push_back(imageDimensions[d]/2);
	}

	float k[3];
	std::vector<int> r(3);
	for(size_t u=0; u<ungriddedData.points(); u++)
	{
    int readoutPoint = u%m_trajectory->numReadoutPoints;
    int readout = u/m_trajectory->numReadoutPoints;
		float densityCompensation;
		trajectoryCoordinates(readoutPoint, readout, m_trajectory, k, &densityCompensation);
		complexFloat ungriddedValue = densityCompensation*ungriddedData.signalValue(u);
		for(size_t g=0; g<griddedData->points(); g++)
		{
			singleToMultiIndex(g, imageDimensions, r);
			float exponentArgument = 0;
			for(int d=0; d<numDimensions(); d++)
			{
				int p = r[d]-imageCenter[d];
				exponentArgument += k[d]*p*axisScale[d];
			}
      exponentArgument *= static_cast<float>(2.0f * M_PI);
			complexFloat griddedValue = complexFloat(std::cos(exponentArgument), std::sin(exponentArgument))*ungriddedValue*imageScale;
			griddedData->addSignalValue(g, griddedValue);
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

std::vector<MRdata *> Gridding::kSpaceToImage(MRdata &ungriddedData, bool returnGriddedKspace)
{
	std::vector<MRdata*> output;
	MRdata* griddedKspace = grid(ungriddedData, Forward);
	MRdata* image;
	if(returnGriddedKspace)
	{
		output.push_back(griddedKspace);
		image = new MRdata(griddedKspace);
	} else {
		image = griddedKspace;
	}
	image->writeToOctave("temp");
	image->fftShift();
	image->fft(FFTW_BACKWARD);
	image->fftShift();

	deapodize(*image);
	image->crop(imageDimensions());

	output.push_back(image);

	return output;
}

MRdata *Gridding::imageToKspace(const MRdata &image)
{
	MRdata griddedKspace = image;
	griddedKspace.pad(m_gridDimensions);
	deapodize(griddedKspace);

	griddedKspace.fftShift();
	griddedKspace.fft(FFTW_FORWARD);
	griddedKspace.fftShift();

	MRdata* ungriddedKspace = grid(griddedKspace, Inverse);

	return ungriddedKspace;
}

int Gridding::numDimensions()
{
  return m_trajectory->numDimensions;
}

std::vector<int> Gridding::gridDimensions()
{
	return m_gridDimensions;
}

