extern "C" 
{
#include "cones.h"
}

#include <mex.h>
#include <math.h>
#include <string.h>
#include <vector>

mxArray* arrayToMxArray(void* array, mxClassID classID, mxComplexity complexity, std::vector<mwSize> size)
{
	int elementSize;
	switch(classID)
	{
		case mxSINGLE_CLASS:
			elementSize = sizeof(float);
			break;
		default:
			exit(1);
	}
	mxArray* mxarray = mxCreateNumericArray(3, size.data(), classID, complexity);
	int points = size[0];
	for(int n=1; n<size.size(); n++)
		points *= size[n];
	memcpy(mxGetData(mxarray), array, points*elementSize); 

	return mxarray;
}

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

	Cones* cones = generateCones(fieldOfViewXY, fieldOfViewZ, NULL, spatialResolutionXY, spatialResolutionZ, 16, 1, NoCompensation, readoutDuration, samplingInterval, filterFieldOfView, maxGradientAmplitude, maxSlewRate, StoreAll);
	Trajectory* trajectory = &cones->trajectory;

	if(!nlhs)
		return;

	const char* infoFields[] = {"fieldOfView", "spatialResolution"};
	mxArray* info = mxCreateStructMatrix(1, 1, 2, infoFields);
	mxSetField(info, 0, "fieldOfView", mxCreateDoubleMatrix(3,1,mxREAL));
	mxSetField(info, 0, "spatialResolution", mxCreateDoubleMatrix(3,1,mxREAL));
	for(int d=0; d<3; d++)
	{
		mxGetPr(mxGetField(info, 0, "fieldOfView"))[d] = cones->trajectory.fieldOfView[d];
		mxGetPr(mxGetField(info, 0, "spatialResolution"))[d] = cones->trajectory.spatialResolution[d];
	}

	const char* outputFields[] = {"info","gradients","kspace"};
	mxArray* output = mxCreateStructMatrix(1, 1, 3, outputFields);
	mxSetField(output, 0, "info", info);

	std::vector<mwSize> gradientDimensions = {trajectory->waveformPoints, 3, trajectory->readouts};
	mxArray* gradientWaveforms = arrayToMxArray(trajectory->gradientWaveforms, mxSINGLE_CLASS, mxREAL, gradientDimensions);
	mxSetField(output, 0, "gradients", gradientWaveforms);

	std::vector<mwSize> trajectoryDimensions = {trajectory->readoutPoints, 3, trajectory->readouts};
	mxArray* kSpaceCoordinates = arrayToMxArray(trajectory->kSpaceCoordinates, mxSINGLE_CLASS, mxREAL, trajectoryDimensions);
	mxSetField(output, 0, "kspace", kSpaceCoordinates);

	plhs[0] = output;

	freeCones(cones);
}

