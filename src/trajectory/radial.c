#include "radial.h"

#include "mrgradient.h"
#include "trajectory.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

struct Trajectory* generateRadial(float* fieldOfView, enum AngleShape thetaShape, float* spatialResolution, enum AngleShape kShape, int dimensions, int fullProjection, float pct, float gmax, float maxSlewRate, float Ts)
{
  struct Trajectory *trajectory = newTrajectory();

	float angleRange;
	float kmax[3];

	if(fullProjection)
		angleRange = M_PI;
	else
		angleRange = 2*M_PI;

  trajectory->numDimensions = dimensions;
  trajectory->numReadoutPoints = 0;
	float spatialResolutionMin = INFINITY;
	int n;
  for(n=0; n<trajectory->numDimensions; n++)
	{
		trajectory->fieldOfView[n] = pct*fieldOfView[n];
		trajectory->spatialResolution[n] = spatialResolution[n];
		adjustSpatialResolution(trajectory->fieldOfView[n], &trajectory->imageDimensions[n], &trajectory->spatialResolution[n]);
		kmax[n] = 10/trajectory->spatialResolution[n];
    trajectory->numReadoutPoints = fmax(ceil(10*trajectory->fieldOfView[n]/trajectory->spatialResolution[n]), trajectory->numReadoutPoints);
		spatialResolutionMin = fmin(trajectory->spatialResolution[n], spatialResolutionMin);
	}

  trajectory->numReadoutPoints += trajectory->numReadoutPoints%2;

	float *angles;
	float *dcfScale;

	float *theta;
	float *phi;

  if(trajectory->numDimensions==2)
    calculateAngles(0, angleRange, thetaShape, trajectory->fieldOfView, kShape, kmax, &angles, NULL, &dcfScale, &trajectory->numReadouts);
	else
    calculateAngles3D(!fullProjection, thetaShape, kShape, trajectory->imageDimensions, &theta, &phi, NULL, &dcfScale, &trajectory->numReadouts);

  for(n=0; n<trajectory->numDimensions; n++)
		kmax[n] *= 2;

//	trajectory->numReadoutPoints = ceil(10*fmax(fovx/resx, fovy/resy));

	trajectory->maxGradientAmplitude = gmax;
	trajectory->maxSlewRate = maxSlewRate;
	trajectory->samplingInterval = Ts;

	float *gradientWaveform;
	int nramp;
  grd_readout(1, 1, trajectory->numReadoutPoints, spatialResolutionMin, gmax, maxSlewRate, Ts, &gradientWaveform, &nramp, &trajectory->numPreReadoutPoints, &trajectory->numWaveformPoints);

  allocateTrajectory(trajectory, trajectory->numReadoutPoints, trajectory->numWaveformPoints, trajectory->numDimensions, 1, trajectory->numReadouts, STORE_ALL);

	int r;
  if(trajectory->numDimensions==2)
	{
		float* nullGradientWaveform = (float*)calloc(trajectory->numWaveformPoints, sizeof(float));
		rotateBasis(gradientWaveform, nullGradientWaveform, trajectory, angleRange);
	}
	else
	{
		float* axisK = (float*)malloc(trajectory->numWaveformPoints*sizeof(float));

		float scale[3];
    for(r=0; r<trajectory->numReadouts; r++)
		{
			scale[0] = sin(theta[r])*cos(phi[r]);
			scale[1] = sin(theta[r])*sin(phi[r]);
			scale[2] = cos(theta[r]);

			int d;
			for(d=0; d<3; d++)
			{
				float *axisGradientWaveform = trajectoryGradientWaveform(trajectory, r, d);
				for(n=0; n<trajectory->numWaveformPoints; n++)
				{
					axisGradientWaveform[n] = scale[d]*gradientWaveform[n];
				}
				gradientToKspace(axisGradientWaveform, axisK, Ts, trajectory->numWaveformPoints);
        memcpy(trajectoryKspaceWaveform(trajectory,r,d), &axisK[trajectory->numPreReadoutPoints], trajectory->numReadoutPoints*sizeof(float));
			}
		}
	}

  for(r=0; r<trajectory->numReadouts; r++)
	{
    for(n=0; n<trajectory->numReadoutPoints; n++)
		{
			float weight = fabs(trajectory->kSpaceCoordinates[n]);
      if(trajectory->numDimensions==3)
				weight *= weight;
      trajectory->densityCompensation[r*trajectory->numReadoutPoints+n] = weight;
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
