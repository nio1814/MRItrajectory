/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#ifndef GRIDDING_H
#define GRIDDING_H

#include "mrdata.h"

#include <vector>

struct Trajectory;

class GriddingTest;

class Gridding
{
	friend class GriddingTest;
public:
  enum Direction{FORWARD, INVERSE};
	Gridding(const Trajectory *trajectory, float oversamplingRatio=1.5, int kernelWidth=4);

	std::vector<int> imageDimensions();

	MRdata* grid(MRdata& inputData, Direction direction);
	MRdata* conjugatePhaseForward(const MRdata& ungriddedData);

	void deapodize(MRdata& oversampledImage);

	std::vector<MRdata*> kSpaceToImage(MRdata& ungriddedData, bool returnGriddedKspace=false);
	MRdata* imageToKspace(const MRdata& image);

private:
	int numDimensions();
	void nearestGriddedPoint(const std::vector<float> &ungriddedPoint, std::vector<float> &griddedPoint);
	float lookupKernelValue(float x);
	std::vector<std::vector<float> > kernelNeighborhood(const std::vector<float> &ungriddedPoint, const std::vector<float> &griddedPoint);

	std::vector<int> gridDimensions();

	const Trajectory* m_trajectory;
	float m_oversamplingFactor;
	int m_kernelWidth;
	std::vector<float> m_kernelLookupTable;
	std::vector<std::vector<float> > m_deapodization;
	std::vector<int> m_gridDimensions;
	std::vector<float> m_normalizedToGridScale;
	float m_coordinateScale;
};

#endif // GRIDDING_H
