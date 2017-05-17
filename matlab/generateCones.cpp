extern "C" 
{
#include "cones.h"
}

#include <mex.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const double *input;

	input = mxGetPr(prhs[0]);
	float fieldOfViewXY = input[0];
	float fieldOfViewZ = input[1];
	input = mxGetPr(prhs[1]);
	float spatialResolutionXY = input[0];
	float spatialResolutionZ = input[1];
	float readoutDuration = *mxGetPr(prhs[2]);

	float filterFieldOfView = fmaxf(fieldOfViewXY, fieldOfViewZ);

	float samplingInterval = 4e-6;
	float maxGradientAmplitude = 4;
	float maxSlewRate = 15000;

	generateCones(fieldOfViewXY, fieldOfViewZ, NULL, spatialResolutionXY, spatialResolutionZ, 16, 1, NoCompensation, readoutDuration, samplingInterval, filterFieldOfView, maxGradientAmplitude, maxSlewRate, StoreAll);
}

