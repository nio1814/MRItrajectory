/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#include "spiral.h"

#include "arrayops.h"
#include "variabledensity.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const float SPIRAL_DESIGN_SAMPLING_INTERVAL = 1e-6;

#ifndef GYROMAGNETIC_RATIO
#define GYROMAGNETIC_RATIO 4257
#endif

float calcc1Spiral(float fieldOfView, float interleaves)
{
	return 4.0f*M_PI*M_PI*fieldOfView*fieldOfView/(interleaves*interleaves);
}

float calcGtwistSpiral(float kr, float fr, float interleaves)
{
    float gtwist;
	float c1 = calcc1Spiral(fr, interleaves);
    /*float c2 = 1.0f;*/	/* only twist when necessary */
    float c2 = 0.0f;	/* always twist */

    gtwist = c1*kr*kr - c2;
    gtwist = sqrt((gtwist>0)?gtwist:0);

    return gtwist;
}

float calcFovFermatFloret(float fieldOfView, float kr, float a0)
{
	return 2*a0*fieldOfView*fieldOfView*kr;
}

int generateSpiral(float fieldOfViewInitial, float spatialResolution, struct VariableDensity *variableDensity, int interleaves, int maxPoints, float samplingInterval, enum SpiralType sptype, float floretAngle, float maxGradientAmplitude, float maxSlewRate, float **kSpaceCoordinates, float **gradientWaveforms, int *readoutPoints)
{
	int oversamplingRatio = round(samplingInterval/SPIRAL_DESIGN_SAMPLING_INTERVAL);
	float gammaDeltaTime = GYROMAGNETIC_RATIO*SPIRAL_DESIGN_SAMPLING_INTERVAL; /* gammaDeltaTime*g = dk */
	float maxGradientDelta = maxSlewRate*SPIRAL_DESIGN_SAMPLING_INTERVAL;	/* the most the gradients can change in 1 raster period */

	float kSpaceMaxRadial = 5/spatialResolution;

	/*float phi = 0;
	float deltaPhi;
    float phiNext;


    float kr;
    float dkr;
    float dkr1, dkr2;
	float deltaPhiDeltaKr;*/
    /*float knorm;*/
	/*int n;
    float *k;
	int *stepSign;
    float gs;
    float gxy;

    float knext[3];
    float dk[3];
    float u[3];
	float projectedGradient;
    float uxym;
	float gradientMagnitudeRange[2];
    float gm;
    float gscale;
	float maxReadoutGradientAmplitude;
    int d;
    float gtwist;
    float term;
	float fov;*/
	int stepBack; /* 1 when backtracking */
	int valid = 1;
    int m;

	int maxDesignPoints = oversamplingRatio*maxPoints;
	float* kDesign = (float*)calloc(2*maxDesignPoints, sizeof(float));
	int* stepSign = (int*)malloc(maxDesignPoints*sizeof(int));

	int n;
	for(n=0; n<maxPoints; n++)
		stepSign[n] = 1;

	kDesign[0] = gammaDeltaTime*maxGradientDelta;
	kDesign[2] = 2*gammaDeltaTime*maxGradientDelta;

	float kr = kDesign[2];

//	dkr1 = kDesign[0];
//	dkr2 = kDesign[2] - kDesign[0];

    n = 1;
    m = 1;

	while((kr <= kSpaceMaxRadial) && valid)
    {
        /*for(d=0; d<ndim; d++)
            km[d] = 1.5*k[n*ndim+d] - 0.5*k[(n-1)*ndim+d];*/
		  float deltaKr = norm2(&kDesign[n*2], 2) - norm2(&kDesign[(n-1)*2], 2);
//		  if(n>=2)
//			dkr2 = norm2(&k[(n-1)*ndim], ndim) - norm2(&k[(n-2)*ndim], ndim);
//		  else
//			  dkr2 = 0;
	/*`	dkr = 2*dkr1 - dkr2;*/
//		dkr = dkr1;

	   /*kmr = norm2(km, ndim);
	   knorm = kr/kSpaceMaxRadial;*/
		  float fieldOfView;
		getFieldOfView(variableDensity, kr, &fieldOfViewInitial, &fieldOfView, 1);
	   float maxReadoutGradientAmplitude = fmin(calculateMaxReadoutGradientAmplitude(fieldOfView, samplingInterval), maxGradientAmplitude);

	   float gtwist;
        if(sptype==spFERMAT)
        {
			fieldOfView = calcFovFermatFloret(fieldOfView, kr, floretAngle);
			gtwist = 2*M_PI*kr*fieldOfView/interleaves;
        }
        else
			gtwist = calcGtwistSpiral(kr, fieldOfView, interleaves);

		float phi = atan2(kDesign[n*2+1], kDesign[n*2]);
		float deltaPhiDeltaKr = gtwist/kr;
        /*dkr = kmr-kr;*/
		float deltaPhi = deltaPhiDeltaKr*deltaKr;
		float phiNext = phi+deltaPhi;

		float kDesignNext[2];
		kDesignNext[0] = cos(phiNext);
		kDesignNext[1] = sin(phiNext);
        /*scalefloats(knext, ndim, kmr);*/
		scalefloats(kDesignNext, 2, kr+deltaKr);

		float deltaK[2];
		subtractArray(kDesignNext, &kDesign[n*2], deltaK, 2);
		float deltaKmagnitude = norm2(deltaK,2);
		float kSpaceDirection[2];
//		copyfloats(deltaK, kSpaceDirection, 2);
//		scalefloats(kSpaceDirection, 2, 1.0f/norm2(deltaK,2));

		float gradient[2];
		int d;
		for(d=0; d<2; d++)
		{
			kSpaceDirection[d] = deltaK[d]/deltaKmagnitude;
			gradient[d] = (kDesign[n*2+d]-kDesign[(n-1)*2+d])/gammaDeltaTime;
		}

//		uxym = sqrt(kSpaceDirection[0]*kSpaceDirection[0]+kSpaceDirection[1]*kSpaceDirection[1]);
		float projectedGradient = dot(kSpaceDirection, gradient, 2);
		float term = (maxGradientDelta*maxGradientDelta - (gradient[0]*gradient[0] + gradient[1]*gradient[1])) + projectedGradient*projectedGradient;

        if(term>=0)
        {
			float s;
			float gradientMagnitudeRange[2];
            for(d=0; d<2; d++)
            {
                if(d==0)
					s = -1.0f;
                else
					s = 1.0f;

				gradientMagnitudeRange[d] = projectedGradient + s*sqrt(term);
            }
	if(gradientMagnitudeRange[0]>gradientMagnitudeRange[1])
				stepBack = 1;
            else
				stepBack = 0;

			if(gradientMagnitudeRange[stepSign[n]]>maxReadoutGradientAmplitude)
				stepBack = 1;
        }
        else
			stepBack = 1;

		if(!stepBack)
        {
			float gradientMagnitude;
			float gradientMagnitudeRange[2];
			if(stepSign[n]==1)
				gradientMagnitude = gradientMagnitudeRange[1];
            else
				gradientMagnitude = gradientMagnitudeRange[0];

			for(d=0; d<2; d++)
            {
				gradient[d] = gradientMagnitude*kSpaceDirection[d];
				kDesign[(n+1)*2+d] = kDesign[n*2+d] + gradient[d]*gammaDeltaTime;
            }
            n++;
        }
        else
        {
			while((n>3) && (stepSign[n-1]==-1))
                n--;
			stepSign[n-1] = -1;
            n -= 2;
        }
		kr = norm2(&kDesign[n*2], 2);

		valid = n<(maxPoints-1);
        m++;
    }

	if(valid)
    {
		*readoutPoints = (n-1)/oversamplingRatio;

		*kSpaceCoordinates = (float*)malloc(*readoutPoints*2*sizeof(float));
		*gradientWaveforms = (float*)malloc(*readoutPoints*2*sizeof(float));

		int d;
		for(d=0; d<2; d++)
        {
			(*kSpaceCoordinates)[*readoutPoints*d] = kDesign[d];
			(*gradientWaveforms)[*readoutPoints*d] = kDesign[d]/gammaDeltaTime;
        }

		for(n=1; n<*readoutPoints; n++)
        {
			for(d=0; d<2; d++)
            {
				(*kSpaceCoordinates)[*readoutPoints*d+n] = kDesign[oversamplingRatio*n*2+d];
				(*gradientWaveforms)[*readoutPoints*d+n] = (kDesign[oversamplingRatio*n*2+d]-kDesign[oversamplingRatio*(n-1)*2+d])/gammaDeltaTime/oversamplingRatio;
            }
        }
    }

	free(kDesign);
	free(stepSign);

	return valid;
}

void generateSpirals(struct VariableDensity *variableDensity, float fieldOfView, float spatialResolution, float readoutDuration, float samplingInterval, int interleavesDesired, enum SpiralType sptype, float floretAngle, float fovFilt, float maxGradientAmplitude, float maxSlewRate)
{
	struct Trajectory *trajectory = (struct Trajectory*)malloc(sizeof(struct Trajectory));
	adjustSpatialResolution(fieldOfView, trajectory->imageDimensions, &spatialResolution);
//    float *g = NULL;
//    float *gbasis;
//    float *gx1, *gx2, *gy1, *gy2, *gz1, *gz2;
//    float *k = NULL;
//    float *kx, *ky;
	float kr;
//    float *grew[3] = {NULL,NULL,NULL};
//    float *g3 = NULL;
//    float fovtemp;
	float kSpaceMaxRadial = 5/spatialResolution;
//	float maxReadoutGradientAmplitude;

//    int Ndes = Tdes/dt;
//    float kstart;
//    int nstart;
	int n;
//    int m;
//    int nrew = 200;

//    int glentemp;
	int readoutPointsTest;
//    int ndim = 2;
//    int ndimOut;
//    int d;

//    float phi;

    /* calculate maximum readout gradient */
    if(fovFilt)
		trajectory->maxReadoutGradientAmplitude = fmin(calculateMaxReadoutGradientAmplitude(fovFilt, samplingInterval), maxGradientAmplitude);
    else
		trajectory->maxReadoutGradientAmplitude = maxGradientAmplitude;

	int readoutPointsDesired = readoutDuration/samplingInterval;
	readoutPointsDesired += readoutPointsDesired%2;

		int interleavesLow = 1;
		int interleavesHigh = 0;

	   for(n=0; n<variableDensity->steps; n++)
	   {
		   /*kr = kSpaceMaxRadial*variableDensity->kn[n];*/
		   float fieldOfViewFinal;
		   kr = variableDensity->step[n].kr;
			getFinalFieldOfView(variableDensity, &fieldOfView, &fieldOfViewFinal, 1);
		   if(sptype==spFERMAT)
			fieldOfViewFinal = calcFovFermatFloret(fieldOfViewFinal, kr, floretAngle);

			interleavesHigh = fmax(interleavesHigh, 2.0f*M_PI*kr*fieldOfViewFinal);
	   }

	   if(kr<kSpaceMaxRadial)
       {
		   kr = kSpaceMaxRadial;
		   float fieldOfViewAtKr;
		   getFinalFieldOfView(variableDensity, &fieldOfView, &fieldOfViewAtKr, 1);
          if(sptype==spFERMAT)
		   fieldOfViewAtKr = calcFovFermatFloret(fieldOfViewAtKr, kr, floretAngle);

		   interleavesHigh = fmax(interleavesHigh, 2.0f*M_PI*kr*fieldOfViewAtKr);
       }

		float* gradientReadoutWaveformsBasis = (float*)calloc(2*readoutPointsDesired, sizeof(float));
		float* gradientWaveformsTest = NULL;
		float* kSpaceCoordinatesTest = NULL;

		while((interleavesHigh-interleavesLow)>1)
        {
			int interleavesTest = (interleavesLow+interleavesHigh)/2.0f;

			if(!readoutDuration)
			{
				interleavesTest = interleavesDesired;
				trajectory->readoutPoints = 1e5;
			}

			if(gradientWaveformsTest)
            {
				free(gradientWaveformsTest);
				gradientWaveformsTest = NULL;
            }

			if(kSpaceCoordinatesTest)
            {
				free(kSpaceCoordinatesTest);
				kSpaceCoordinatesTest = NULL;
            }

		  if(generateSpiral(fieldOfView, spatialResolution, variableDensity, interleavesTest, trajectory->readoutPoints, samplingInterval, sptype, floretAngle, trajectory->maxReadoutGradientAmplitude, maxSlewRate, &kSpaceCoordinatesTest, &gradientWaveformsTest, &readoutPointsTest))
            {

				if(readoutPointsTest>trajectory->readoutPoints)
					interleavesLow = interleavesTest;
                else
                {
					interleavesHigh = interleavesTest;
					trajectory->readoutPoints = readoutPointsTest;
					memcpy(gradientReadoutWaveformsBasis, gradientWaveformsTest, 2*trajectory->readoutPoints*sizeof(float));
					trajectory->readouts = interleavesTest;
                }
            }
            else
				interleavesLow = interleavesTest;

			printf(" %d interleaves %d pts\n", interleavesTest, readoutPointsTest);

			if(!readoutDuration)
				interleavesHigh = interleavesLow;
        }

	float* gradientWaveformsBasis[2];
	traverseKspaceToZero(gradientReadoutWaveformsBasis, &gradientReadoutWaveformsBasis[trajectory->readoutPoints], NULL, trajectory->readoutPoints, samplingInterval, maxGradientAmplitude, maxSlewRate, &gradientWaveformsBasis[0], &gradientWaveformsBasis[1], NULL, &trajectory->waveformPoints);

	allocateTrajectory(trajectory, trajectory->readoutPoints, trajectory->waveformPoints, 2, 1, trajectory->readouts, StoreAll);
	trajectory->imageDimensions[1] = trajectory->imageDimensions[0];
	for(n=0; n<2; n++)
	{
		trajectory->spatialResolution[n] = spatialResolution;
		trajectory->fieldOfView[n] = fieldOfView;
	}
	trajectory->samplingInterval = samplingInterval;
	trajectory->variableDensity = variableDensity;
	trajectory->maxGradientAmplitude = maxGradientAmplitude;
	trajectory->maxSlewRate = maxSlewRate;

	for(n=0; n<trajectory->readouts; n++)
    {
		float phi = (2.0f*M_PI*n)/trajectory->readouts;

		float* gx = &trajectory->gradientWaveforms[2*n*trajectory->waveformPoints];
		float* gy = &trajectory->gradientWaveforms[(2*n+1)*trajectory->waveformPoints];

		float* kx = &trajectory->kSpaceCoordinates[2*n*trajectory->readoutPoints];
		float* ky = &trajectory->kSpaceCoordinates[(2*n+1)*trajectory->readoutPoints];

		memcpy(gx, gradientWaveformsBasis[0], trajectory->waveformPoints*sizeof(float));
		memcpy(gy, gradientWaveformsBasis[1], trajectory->waveformPoints*sizeof(float));
		scalecomplex(gx, gy, cos(phi), sin(phi), trajectory->waveformPoints);

		gradientToKspace(gx, kx, samplingInterval, trajectory->readoutPoints);
		gradientToKspace(gy, ky, samplingInterval, trajectory->readoutPoints);

		calcSpiralDcf(gx, gy, kx, ky, trajectory->readoutPoints, &trajectory->densityCompensation[n*trajectory->readoutPoints]);
    }

    return;
}


