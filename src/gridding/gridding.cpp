/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#include "gridding.h"

#include "trajectory.h"

extern "C" {
#include "numericalrecipes.h"
#include "arrayops.h"
}

Gridding::Gridding(const Trajectory *trajectory) :
	m_trajectory(trajectory), m_oversamplingFactor(1.5)
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

	std::vector<int> gridDims = gridDimensions();
	for(int d=0; d<m_trajectory->dimensions; d++)
	{
		for(int n=0; n<gridDims[d]; n++)
		{
			float arg = M_PI*m_kernelWidth*(-gridDims[d]/2.0f + n)/gridDims[d];
			arg *= arg;
			arg -= beta*beta;

			float value;
			if(arg<0)
				value = -1.0f/(2.0f*sqrtf(-arg))*(expf(-sqrtf(-arg))-expf(sqrtf(-arg)));
			else
				value = sinf(sqrtf(arg))/sqrtf(arg);
			value = scale/value;
			m_deapodization.push_back(value);
		}
	}
}

MRdata *Gridding::grid(const MRdata &ungriddedData)
{

}

std::vector<int> Gridding::gridDimensions()
{
	std::vector<int> dims;

	for(int d=0; d<m_trajectory->dimensions; d++)
	{
		int dimension = m_oversamplingFactor*m_trajectory->imageDimensions[d];
		dimension += dimension%2;
		dims.push_back(dimension);
	}

	return dims;
}


