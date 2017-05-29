#include "trajectory.h"

#include "arrayops.h"
#include "mathops.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

const float GYROMAGNETIC_RATIO = 4257.59;
const float REWIND_SAMPLING_INTERVAL = 1e-6;
const int MAX_REWIND_POINTS = 4096;

void initializeTrajectory(struct Trajectory *trajectory)
{
	trajectory->densityCompensation = NULL;
	trajectory->kSpaceCoordinates = NULL;
	trajectory->gradientWaveforms = NULL;
	trajectory->gradientWaveformsShort = NULL;
}

void adjustSpatialResolution(float fieldOfView, int *imageDimension, float *spatialResolution)
{
	*imageDimension = ceil(10.0*fieldOfView/(*spatialResolution));
	*imageDimension += (*imageDimension)%2;

	*spatialResolution = 10.0f*fieldOfView/ *imageDimension;
}

void gradientToKspace(float* gradient, float* kSpaceCoordinate, float samplingInterval, int length)
{
    int i;
    float gradientSum = 0.0f;	/* accumulation of gradient */

    for(i=0; i<length; i++)
    {
        gradientSum += gradient[i]*GYROMAGNETIC_RATIO*samplingInterval;
        kSpaceCoordinate[i] = gradientSum;
    }

    return;
}

void gradientWaveformToShort(const float *gradientsFloat, int points, float maxGradientAmplitudeScanner, short* gradientsShort)
{
	float gaussPerCMTo16bit = 32766/maxGradientAmplitudeScanner;

	int n=0;

	for(n=0; n<points-1; n++)
		gradientsShort[n] = (short)(gradientsFloat[n]*gaussPerCMTo16bit) & 0xfffe;

	gradientsShort[n] = (short)(gradientsFloat[n]*gaussPerCMTo16bit) | 0x0001;
}

void allocateTrajectory(struct Trajectory *trajectory, int readoutPoints, int waveformPoints, int dimensions, int bases, int readouts, enum WaveformStorageType storage)
{
	int count;

	if(storage==StoreBasis)
		count = bases;
	else
		count = readouts;

	if(waveformPoints)
		trajectory->gradientWaveforms = (float*)malloc(dimensions*count*waveformPoints*sizeof(float));
	if(readoutPoints)
	{
		trajectory->kSpaceCoordinates = (float*)malloc(dimensions*count*readoutPoints*sizeof(float));
		trajectory->densityCompensation = (float*)malloc(count*readoutPoints*sizeof(float));
	}

	trajectory->storage = storage;
	trajectory->readouts = readouts;
	trajectory->bases = bases;
	trajectory->readoutPoints = readoutPoints;
	trajectory->waveformPoints = waveformPoints;
	trajectory->dimensions = dimensions;
}

void traverseKspace(float *gradientInitial, float *kSpaceCoordinatesInitial, int dimensions, float samplingInterval, float *kSpaceCoordinatesFinal, float maxGradientAmplitude, float maxSlewRate, float **gradientRewindX, float **gradientRewindY, float **gradientRewindZ, int *pointsRewind)
{
	int oversamplingRatio = roundToInteger(samplingInterval/REWIND_SAMPLING_INTERVAL);
	float gain = 20;

	 /*float kacc = .0001;		 Accuracy with which to rewind k.
	 float gacc = 4.0/50000.0;	 Accuracy with which to rewind gradient. */

/* 	Set precission so final point can be reach with one additional step */
	float maxGradientDelta = maxSlewRate*REWIND_SAMPLING_INTERVAL;
	float gacc = oversamplingRatio*maxGradientDelta;
	float kacc = oversamplingRatio*gacc*maxGradientDelta;

	float *gradient = (float*)malloc(MAX_REWIND_POINTS*3*sizeof(float));
	float *karr = (float*)malloc(MAX_REWIND_POINTS*sizeof(float));
	float kSpaceCoordinates[3], kret[3], gradientNormalized[3], gtarg[3], deltaGradientvect[3], deltaGradient[3];
	int keepgoing = 1;
	int n;
	int copySize = dimensions*sizeof(float);
	float a;

	/* First rotate trajectory. This is to try to get design to be insensitive to starting rotation.*/
	memcpy(kSpaceCoordinates, kSpaceCoordinatesInitial, copySize);
	memcpy(gradient, gradientInitial, copySize);

	n = 0;

	memcpy(kret, kSpaceCoordinates, copySize);

	while(keepgoing==1 && n<(MAX_REWIND_POINTS-1))
	{
		/*kret = k + GYROMAGNETIC_RATIO*0.5*g[n]/maxSlewRate;*/
		memcpy(kret, &gradient[dimensions*n], copySize);
		scalefloats(kret, dimensions, GYROMAGNETIC_RATIO*0.5f*norm2(&gradient[dimensions*n],dimensions)/maxSlewRate);
		addfloats(kret, kSpaceCoordinates, kret, dimensions);

		if (distance(kret,kSpaceCoordinatesFinal,dimensions) < kacc)	/*% Then keep ramping toward zero. */
		{
		   /*within k accuracy*/
			/*gradientNormalized = g[n];*/
			memcpy(gradientNormalized, &gradient[dimensions*n], copySize);
			scalefloats(gradientNormalized, dimensions, 1/norm2(gradientNormalized,dimensions));

			if (norm2(&gradient[dimensions*n],dimensions) > maxGradientDelta)
			{
				/*g[n+1] = g[n] - gradientNormalized * maxSlewRate * samplingInterval;*/
				scalefloats(gradientNormalized, dimensions, maxGradientDelta);
				subtractArray(&gradient[dimensions*n], gradientNormalized, &gradient[dimensions*(n+1)], dimensions);
			}
			else
				/*g[n+1] = 0;*/
				memset(&gradient[dimensions*(n+1)], 0, copySize);

			/*gtarg = 0;*/
			memset(gtarg, 0, copySize);
		}
		else
		{
			/*gtarg = gain*(kSpaceCoordinatesFinal-kret);*/
			subtractArray(kSpaceCoordinatesFinal, kret, gtarg, dimensions);
			scalefloats(gtarg, dimensions, gain);

			a = norm2(gtarg, dimensions);
			if (a > maxGradientAmplitude)
				/*gtarg*maxGradientAmplitude/norm2(gtarg, dimensions);*/
				scalefloats(gtarg, dimensions, maxGradientAmplitude/a);

			/*deltaGradientvect = gtarg - g[n];*/
			subtractArray(gtarg, &gradient[dimensions*n], deltaGradientvect, dimensions);

			a = norm2(deltaGradientvect,dimensions);
			if (a > (maxGradientDelta))
			{
				/*deltaGradient = deltaGradientvect*maxSlewRate*samplingInterval;*/
				memcpy(deltaGradient, deltaGradientvect, copySize)   ;
				scalefloats(deltaGradient, dimensions, maxGradientDelta/a);
			}
			else
			 /*   deltaGradient = deltaGradientvect;*/
				memcpy(deltaGradient, deltaGradientvect, copySize);

			/*g[n+1] = g[n] + deltaGradient;*/
			addfloats(&gradient[dimensions*n], deltaGradient, &gradient[dimensions*(n+1)], dimensions);
		}

		karr[n] = distance(kSpaceCoordinatesFinal, kSpaceCoordinates, dimensions);

		memcpy(deltaGradient, &gradient[dimensions*(n+1)], dimensions*sizeof(float));
		scalefloats(deltaGradient, dimensions, GYROMAGNETIC_RATIO*REWIND_SAMPLING_INTERVAL);

		addfloats(kSpaceCoordinates, deltaGradient, kSpaceCoordinates, dimensions);

		n++;

		if(n > MAX_REWIND_POINTS)
	   {
			keepgoing = 0;
		  fprintf(stderr, "goTok1: Error, rewinder length exceeded\n");
	   }
		else if ((norm2(kSpaceCoordinates,dimensions) < kacc) & (norm2(&gradient[dimensions*n], dimensions) < gacc) && !(n%oversamplingRatio))
		/*else if ((norm2(k,dimensions) < kacc) & (norm2(g[n], dimensions) < gacc))*/
			keepgoing = 0;
	}

//	ndes = n;
	*pointsRewind = (n+1)/oversamplingRatio;

	copySize = *pointsRewind*sizeof(float);

	*gradientRewindX = (float*)malloc(copySize);
	if(gradientRewindY)
		*gradientRewindY = (float*)malloc(copySize);
	if(gradientRewindZ)
		*gradientRewindZ = (float*)malloc(copySize);

	for(n=1; n<=*pointsRewind; n++)
	{
		(*gradientRewindX)[n-1] = gradient[3*oversamplingRatio*n];
		if(gradientRewindY)
			(*gradientRewindY)[n-1] = gradient[3*oversamplingRatio*n+1];
		if(gradientRewindZ)
			(*gradientRewindZ)[n-1] = gradient[3*oversamplingRatio*n+2];
	}

	n = *pointsRewind-1;

	if(gradientRewindX)
	{
		a = (*gradientRewindX)[n]/samplingInterval;
		if(fabs(a)>maxSlewRate)
		   fprintf(stderr, "gx %f out of range", (*gradientRewindX)[n]);
	}

	if(gradientRewindY)
	{
		a = (*gradientRewindY)[n]/samplingInterval;
		if(fabs(a)>maxSlewRate)
		   fprintf(stderr, "gy %f out of range", (*gradientRewindY)[n]);
	}

	if(gradientRewindZ)
	{
		a = (*gradientRewindZ)[n]/samplingInterval;
		if(fabs(a)>maxSlewRate)
		   fprintf(stderr, "gz %f out of range", (*gradientRewindZ)[n]);
	}

	return;
}

void traverseKspaceFromWaveform(float *gradientOriginalX, float *gradientOriginalY, float *gradientOriginalZ, int pointsOriginal, float *kSpaceCoordinatesFinal, float samplingInterval, float maxGradientAmplitude, float maxSlewRate, float **gradientRewoundX, float**gradientRewoundY, float **gradientRewoundZ, int *pointsRewound)
{
	int dimensions;
	float kSpaceCoordinatesInitial[3] = {0,0,0};
	float gradientInitial[3] = {0,0,0};
	float *gradientRewindX=NULL, *gradientRewindY=NULL, *gradientRewindZ=NULL;
	float **gxp, **gyp, **gzp;
	int pointsRewind;

	if(gradientOriginalX)
	{
		gradientInitial[0] = gradientOriginalX[pointsOriginal-1];
		kSpaceCoordinatesInitial[0] = sumfloats(gradientOriginalX, pointsOriginal)*GYROMAGNETIC_RATIO*samplingInterval;
		dimensions = 1;
	   gxp = &gradientRewindX;
	}
	else
		gxp = NULL;

	if(gradientOriginalY)
	{
		gradientInitial[1] = gradientOriginalY[pointsOriginal-1];
		kSpaceCoordinatesInitial[1] = sumfloats(gradientOriginalY, pointsOriginal)*GYROMAGNETIC_RATIO*samplingInterval;
		dimensions = 2;
	   gyp = &gradientRewindY;
	}
	else
		gyp = NULL;

	if(gradientOriginalZ)
	{
		gradientInitial[2] = gradientOriginalZ[pointsOriginal-1];
		kSpaceCoordinatesInitial[2] = sumfloats(gradientOriginalZ, pointsOriginal)*GYROMAGNETIC_RATIO*samplingInterval;
		dimensions = 3;
	   gzp = &gradientRewindZ;
	}
	else
		gzp = NULL;

	traverseKspace(gradientInitial, kSpaceCoordinatesInitial, dimensions, samplingInterval, kSpaceCoordinatesFinal, maxGradientAmplitude, maxSlewRate, gxp, gyp, gzp, &pointsRewind);

	*pointsRewound = pointsOriginal+pointsRewind;

	if(gradientOriginalX)
	{
		*gradientRewoundX = (float*)malloc(*pointsRewound*sizeof(float));
		memcpy(*gradientRewoundX, gradientOriginalX, pointsOriginal*sizeof(float));
		memcpy(&(*gradientRewoundX)[pointsOriginal], gradientRewindX, pointsRewind*sizeof(float));
	}

	if(gradientOriginalY)
	{
		*gradientRewoundY = (float*)malloc(*pointsRewound*sizeof(float));
		memcpy(*gradientRewoundY, gradientOriginalY, pointsOriginal*sizeof(float));
		memcpy(&(*gradientRewoundY)[pointsOriginal], gradientRewindY, pointsRewind*sizeof(float));
	}

	if(gradientOriginalZ)
	{
		*gradientRewoundZ = (float*)malloc(*pointsRewound*sizeof(float));
		memcpy(*gradientRewoundZ, gradientOriginalZ, pointsOriginal*sizeof(float));
		memcpy(&(*gradientRewoundZ)[pointsOriginal], gradientRewindZ, pointsRewind*sizeof(float));
	}

	free(gradientRewindX);
	free(gradientRewindY);
	free(gradientRewindZ);

	return;
}

void traverseKspaceToZero(float *gradientOriginalX, float *gradientOriginalY, float *gradientOriginalZ, int pointsOriginal, float samplingInterval, float maxGradientAmplitude, float maxSlewRate, float **gradientRewoundX, float**gradientRewoundY, float **gradientRewoundZ, int *pointsRewound)
{
	float kSpaceCoordinatesFinal[3] = {0,0,0};

	return traverseKspaceFromWaveform(gradientOriginalX, gradientOriginalY, gradientOriginalZ, pointsOriginal, kSpaceCoordinatesFinal, samplingInterval, maxGradientAmplitude, maxSlewRate, gradientRewoundX, gradientRewoundY, gradientRewoundZ, pointsRewound);
}

#define SAVE_GRADIENT_DESCRIPTION_LENGTH 256
int saveGradientWaveforms(const char *filename, const float* grad, short dimensions, short interleaves, short points, int readoutPoints, float FOV, float maxGradientAmplitude, float maxGradientAmplitudeScanner, float samplingInterval, const char* description, enum Endian endian)
{
	int status = 1;
	FILE* file;
	double params[11];
	short np;
	short* nw;
	char descriptionSave[SAVE_GRADIENT_DESCRIPTION_LENGTH];
	int d, i;
	short* waveformsI;
	int swapEndian = needEndianSwap(endian);

	file = fopen(filename, "wb");
	if(!file)
	   fprintf(stderr, "Error opening %s for read", filename);
	else
	{
		nw = (short*)malloc(dimensions*sizeof(short));
		for(d=0; d<dimensions; d++)
			nw[d] = interleaves;
		params[0] = 0;	/*unsure what this is for*/
		params[1] = FOV;
		params[2] = interleaves*dimensions;
		params[3] = maxGradientAmplitude;
		params[4] = points;
		params[5] = samplingInterval*1e6;
		params[6] = readoutPoints;
		params[7] = samplingInterval*1e6;
		params[8] = 0;
		params[9] = 0;
		params[10] = 0;
		np = 11;

		memset(descriptionSave, 0, SAVE_GRADIENT_DESCRIPTION_LENGTH*sizeof(char));
		memcpy(descriptionSave, description, strlen(description)*sizeof(char));

		/* Convert G to 16 bit even numbers */
		waveformsI = (short*)malloc(dimensions*interleaves*points*sizeof(short));
		for(d=0; d<dimensions; d++)
			for(i=0; i<interleaves; i++)
				gradientWaveformToShort(&grad[(i*dimensions+d)*points], points, maxGradientAmplitudeScanner, &waveformsI[(d*interleaves+i)*points]);

		fwrite(descriptionSave, sizeof(char), SAVE_GRADIENT_DESCRIPTION_LENGTH, file);

		if(swapEndian)
		{
			swapArrayEndian(&points, 1, sizeof(short));
			swapArrayEndian(&dimensions, 1, sizeof(short));
		}

		fwrite(&points, sizeof(short), 1, file);
		fwrite(&dimensions, sizeof(short), 1, file);

		if(swapEndian)
		{
			swapArrayEndian(&dimensions, 1, sizeof(short));
			swapArrayEndian(nw, dimensions, sizeof(short));
			swapArrayEndian(&np, 1, sizeof(short));
		}

		fwrite(nw, sizeof(short), dimensions, file);
		fwrite(&np, sizeof(short), 1, file);

		if(swapEndian)
		{
			swapArrayEndian(&np, 1, sizeof(short));
			swapArrayEndian(params, np, sizeof(double));
			swapArrayEndian(&points, 1, sizeof(short));
		}

		fwrite(params, sizeof(double), np, file);

		if(swapEndian)
			swapArrayEndian(waveformsI, dimensions*interleaves*points, sizeof(short));
		fwrite(waveformsI, sizeof(short), dimensions*interleaves*points, file);

		fclose(file);
		free(nw);
		free(waveformsI);

		status = 0;
	}

	return status;
}

int numTrajectoryWaveforms(const struct Trajectory *trajectory)
{
	return trajectory->storage==StoreBasis ? trajectory->bases : trajectory->readouts;
}

void trajectoryCoordinates(int readoutPoint, int readout, const struct Trajectory *trajectory, float* coordinates)
{
	int d;
	for(d=0; d<trajectory->dimensions; d++)
		coordinates[d] = trajectory->kSpaceCoordinates[(trajectory->dimensions*readout+d)*trajectory->readoutPoints+readoutPoint];
}

void setTrajectoryPoint(int readoutPoint, int readout, struct Trajectory *trajectory, const float *coordinates, float densityCompensation)
{
	int d;
	for(d=0; d<trajectory->dimensions; d++)
		trajectory->kSpaceCoordinates[(trajectory->dimensions*readout+d)*trajectory->readoutPoints+readoutPoint] = coordinates[d];
	trajectory->densityCompensation[readout*trajectory->readoutPoints+readoutPoint] = densityCompensation;
}
