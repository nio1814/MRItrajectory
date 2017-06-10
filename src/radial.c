/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#include "radial.h"

#include "mrgradient.h"
#include "trajectory.h"
#include "densitycompensation2D.h"

#include <math.h>
#include <stdlib.h>

struct Trajectory* generateRadial(float fovx, float fovy, enum AngleShape thetaShape, float resx, float resy, enum AngleShape kShape, int fullProjection, float pct, float gmax, float maxSlewRate, float Ts)
{
	struct Trajectory *trajectory = (struct Trajectory*)malloc(sizeof(struct Trajectory));
	initializeTrajectory(trajectory);

	float angleRange;
	float kmax[2];

	if(fullProjection)
		angleRange = M_PI;
	else
		angleRange = 2*M_PI;

	trajectory->fieldOfView[0] = pct*fovx;
	trajectory->fieldOfView[1] = pct*fovy;

	trajectory->spatialResolution[0] = resx;
	trajectory->spatialResolution[1] = resy;

	int n;
	for(n=0; n<2; n++)
	{
		adjustSpatialResolution(trajectory->fieldOfView[n], &trajectory->imageDimensions[n], &trajectory->spatialResolution[n]);
		kmax[n] = 10/trajectory->spatialResolution[n];
	}

	float *angles;
	float *dcfScale;
	calculateAngles(0, angleRange, thetaShape, trajectory->fieldOfView, kShape, kmax, &angles, NULL, &dcfScale, &trajectory->readouts);

	for(n=0; n<2; n++)
		kmax[n] *= 2;

	trajectory->readoutPoints = ceil(10*fmax(fovx/resx, fovy/resy));
	trajectory->readoutPoints += trajectory->readoutPoints%2;
	float spatialResolutionMin = fmin(resx, resy);

	trajectory->maxGradientAmplitude = gmax;
	trajectory->maxSlewRate = maxSlewRate;
	trajectory->samplingInterval = Ts;

	float *gradientWaveform;
//	int ndep;
	int nramp;
	grd_readout(1, 1, trajectory->readoutPoints, spatialResolutionMin, gmax, maxSlewRate, Ts, &gradientWaveform, &nramp, &trajectory->preReadoutPoints, &trajectory->waveformPoints);

	allocateTrajectory(trajectory, trajectory->readoutPoints, trajectory->waveformPoints, 2, 1, trajectory->readouts, StoreAll);

	float* nullGradientWaveform = (float*)calloc(trajectory->waveformPoints, sizeof(float));
	rotateBasis(gradientWaveform, nullGradientWaveform, trajectory, angleRange);

	int r;
	for(r=0; r<trajectory->readouts; r++)
	{
		for(n=0; n<trajectory->readoutPoints; n++)
			trajectory->densityCompensation[r*trajectory->readoutPoints+n] = fabs(trajectory->kSpaceCoordinates[n]);
	}

	convolutionDensityCompensation(trajectory, 0, NULL, 2, 2, trajectory->densityCompensation);

	return trajectory;
}
