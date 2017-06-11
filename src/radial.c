#include "radial.h"

#include "mrgradient.h"
#include "trajectory.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

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

	float *theta;
	float *phi;

	if(trajectory->dimensions==2)
		calculateAngles(0, angleRange, thetaShape, trajectory->fieldOfView, kShape, kmax, &angles, NULL, &dcfScale, &trajectory->readouts);
	else
		calculateAngles3D(!fullProjection, thetaShape, kShape, trajectory->imageDimensions, &theta, &phi, NULL, &dcfScale, &trajectory->readouts);

	for(n=0; n<trajectory->dimensions; n++)
		kmax[n] *= 2;

//	trajectory->readoutPoints = ceil(10*fmax(fovx/resx, fovy/resy));

	trajectory->maxGradientAmplitude = gmax;
	trajectory->maxSlewRate = maxSlewRate;
	trajectory->samplingInterval = Ts;

	float *gradientWaveform;
	int nramp;
	grd_readout(1, 1, trajectory->readoutPoints, spatialResolutionMin, gmax, maxSlewRate, Ts, &gradientWaveform, &nramp, &trajectory->preReadoutPoints, &trajectory->waveformPoints);

	allocateTrajectory(trajectory, trajectory->readoutPoints, trajectory->waveformPoints, trajectory->dimensions, 1, trajectory->readouts, StoreAll);

	int r;
	if(trajectory->dimensions==2)
	{
		float* nullGradientWaveform = (float*)calloc(trajectory->waveformPoints, sizeof(float));
		rotateBasis(gradientWaveform, nullGradientWaveform, trajectory, angleRange);
	}
	else
	{
		float* axisK = (float*)malloc(trajectory->waveformPoints*sizeof(float));

		float scale[3];
		for(r=0; r<trajectory->readouts; r++)
		{
			scale[0] = sin(theta[r])*cos(phi[r]);
			scale[1] = sin(theta[r])*sin(phi[r]);
			scale[2] = cos(theta[r]);

			int d;
			for(d=0; d<3; d++)
			{
				float *axisGradientWaveform = trajectoryGradientWaveform(trajectory, r, d);
				for(n=0; n<trajectory->waveformPoints; n++)
				{
					axisGradientWaveform[n] = scale[d]*gradientWaveform[n];
				}
				gradientToKspace(axisGradientWaveform, axisK, Ts, trajectory->waveformPoints);
				memcpy(trajectoryKspaceWaveform(trajectory,r,d), &axisK[trajectory->preReadoutPoints], trajectory->readoutPoints*sizeof(float));
			}
		}
	}

	for(r=0; r<trajectory->readouts; r++)
	{
		for(n=0; n<trajectory->readoutPoints; n++)
		{
			float weight = fabs(trajectory->kSpaceCoordinates[n]);
			if(trajectory->dimensions==3)
				weight *= weight;
			trajectory->densityCompensation[r*trajectory->readoutPoints+n] = weight;
		}

	}

	return trajectory;
}

struct Trajectory* generateRadial2D(float fovx, float fovy, enum AngleShape thetaShape, float resx, float resy, enum AngleShape kShape, int fullProjection, float pct, float gmax, float maxSlewRate, float Ts)
{
	float fieldOfView[3] = {fovx, fovy};
	float spatialResolution[3] = {resx, resy};

	return generateRadial(fieldOfView, thetaShape, spatialResolution, kShape, 2, fullProjection, pct, gmax, maxSlewRate, Ts);
}

struct Trajectory *generateRadial3D(float fovx, float fovy, float fieldOfViewZ, enum AngleShape thetaShape, enum AngleShape phiShape, float resx, float resy, float spatialResolutionZ, int fullProjection, float pct, float gmax, float maxSlewRate, float Ts)
{
	float fieldOfView[3] = {fovx, fovy, fieldOfViewZ};
	float spatialResolution[3] = {resx, resy, spatialResolutionZ};

	return generateRadial(fieldOfView, thetaShape, spatialResolution, phiShape, 3, fullProjection, pct, gmax, maxSlewRate, Ts);
}
