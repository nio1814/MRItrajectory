#include "radial.h"

#include "mrgradient.h"
#include "trajectory.h"

#include <math.h>
#include <stdlib.h>

struct Trajectory* generateRadial(float* fieldOfView, enum AngleShape thetaShape, float* spatialResolution, enum AngleShape kShape, int dimensions, int fullProjection, float pct, float gmax, float maxSlewRate, float Ts)
{
	struct Trajectory *trajectory = (struct Trajectory*)malloc(sizeof(struct Trajectory));
	initializeTrajectory(trajectory);

	float angleRange;
	float kmax[3];

	if(fullProjection)
		angleRange = M_PI;
	else
		angleRange = 2*M_PI;

	trajectory->dimensions = dimensions;
	trajectory->readoutPoints = 0;
	float spatialResolutionMin = INFINITY;
	int n;
	for(n=0; n<trajectory->dimensions; n++)
	{
		trajectory->fieldOfView[n] = pct*fieldOfView[n];
		trajectory->spatialResolution[n] = spatialResolution[n];
		adjustSpatialResolution(trajectory->fieldOfView[n], &trajectory->imageDimensions[n], &trajectory->spatialResolution[n]);
		kmax[n] = 10/trajectory->spatialResolution[n];
		trajectory->readoutPoints = fmax(ceil(10*trajectory->fieldOfView[n]/trajectory->spatialResolution[n]), trajectory->readoutPoints);
		spatialResolutionMin = fmin(trajectory->spatialResolution[n], spatialResolutionMin);
	}

	trajectory->readoutPoints += trajectory->readoutPoints%2;

	float *angles;
	float *dcfScale;

	if(trajectory->dimensions==2)
		calculateAngles(0, angleRange, thetaShape, trajectory->fieldOfView, kShape, kmax, &angles, NULL, &dcfScale, &trajectory->readouts);
	else


	for(n=0; n<2; n++)
		kmax[n] *= 2;

//	trajectory->readoutPoints = ceil(10*fmax(fovx/resx, fovy/resy));

	trajectory->maxGradientAmplitude = gmax;
	trajectory->maxSlewRate = maxSlewRate;
	trajectory->samplingInterval = Ts;

	float *gradientWaveform;
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

	return trajectory;
}

struct Trajectory* generateRadial2D(float fovx, float fovy, enum AngleShape thetaShape, float resx, float resy, enum AngleShape kShape, int fullProjection, float pct, float gmax, float maxSlewRate, float Ts)
{
	float fieldOfView[3] = {fovx, fovy};
	float spatialResolution[3] = {resx, resy};

	return generateRadial(fieldOfView, thetaShape, spatialResolution, kShape, 2, fullProjection, pct, gmax, maxSlewRate, Ts);
}
