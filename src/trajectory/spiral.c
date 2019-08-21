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
#include "mathops.h"

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

	int stepBack; 
	int valid = 1;

	int maxPointsDesign = oversamplingRatio*maxPoints;
	float* kDesign = (float*)calloc(2*maxPointsDesign, sizeof(float));
	int* stepSign = (int*)malloc(maxPointsDesign*sizeof(int));

	int n;
	for(n=0; n<maxPointsDesign; n++)
		stepSign[n] = 1;

	kDesign[0] = gammaDeltaTime*maxGradientDelta;
	kDesign[2] = 2*gammaDeltaTime*maxGradientDelta;

	float kr = kDesign[2];

    n = 1;

	while(kr <= kSpaceMaxRadial && n<(maxPointsDesign-1))
    {
		  float deltaKr = norm2(&kDesign[n*2], 2) - norm2(&kDesign[(n-1)*2], 2);
		  float fieldOfView = fieldOfViewInitial;
		  if(variableDensity)
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
		float deltaPhi = deltaPhiDeltaKr*deltaKr;
		float phiNext = phi+deltaPhi;

		float kDesignNext[2];
		kDesignNext[0] = cos(phiNext);
		kDesignNext[1] = sin(phiNext);
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

		float gradientMagnitudeRange[2];
        if(term>=0)
        {
			float s;

            for(d=0; d<2; d++)
            {
                if(d==0)
					s = -1.0f;
                else
					s = 1.0f;

				gradientMagnitudeRange[d] = projectedGradient + s*sqrt(term);
            }
			gradientMagnitudeRange[0] = fmax(gradientMagnitudeRange[0], 0);
		if(gradientMagnitudeRange[0]>gradientMagnitudeRange[1])
				stepBack = 1;
            else
				stepBack = 0;

		int stepIndex = (stepSign[n]+1)/2;
			if(gradientMagnitudeRange[stepIndex]>.999*maxReadoutGradientAmplitude)
				stepBack = 1;
        }
        else
			stepBack = 1;

		if(!stepBack)
        {
			float gradientMagnitude;
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
    }

	valid = n<(maxPointsDesign-1);

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

	return !valid;
}

void calcSpiralDcf(float *gx, float *gy, float *kx, float *ky, int rolen, float *denscomp)
{
	float *kr, *gradientMagnitude;
	float* Gtwist;
	float *denscomp_it, *denscomp_is;
	float *kPhasexy, *kMagxy;
	int n;

	/* Calculate k-space and gradient magnitudes */
	kr = (float*)malloc(rolen*sizeof(float));
	gradientMagnitude = (float*)malloc(rolen*sizeof(float));
	for(n=0; n<rolen; n++)
	{
        kr[n] = sqrt(kx[n]*kx[n] + ky[n]*ky[n]);
        gradientMagnitude[n] = sqrt(gx[n]*gx[n] + gy[n]*gy[n]);
	}

	Gtwist = (float*)malloc(rolen*sizeof(float));
	kPhasexy = (float*)malloc(rolen*sizeof(float));
	kMagxy = (float*)malloc(rolen*sizeof(float));

	calculatePhase(kx, ky, kPhasexy, rolen, 0, 1.0);
	calculateMagnitude(kx, ky, kMagxy, rolen);
	unwrapPhase(kPhasexy, Gtwist, rolen, rolen, M_PI);
	for(n=0; n<rolen-1; n++)
        Gtwist[n] = fmax((Gtwist[n+1]-Gtwist[n])/(kr[n+1]-kr[n])*kr[n], 0.0f);
	
	Gtwist[rolen-1] = Gtwist[rolen-2];

	/*% Density compensation due to inter-trajectory spacing (ignoring NINT)
	denscomp_it = abs(kcxy)./sqrt(1+Gtwist.^2);*/
	denscomp_it = (float*)malloc(rolen*sizeof(float));
	for(n=0; n<rolen; n++)
    {
        if(Gtwist[n]>=0)
            denscomp_it[n] = kMagxy[n]/sqrt(1+Gtwist[n]*Gtwist[n]);
        else
            denscomp_it[n] = 0;
    }

	/*% Density compensation due to inter-sample spacing
	denscomp_is = (gr+[gr(2:end); gr(end)])/2;*/
	denscomp_is = (float*)malloc(rolen*sizeof(float));
	addfloats(gradientMagnitude, &(gradientMagnitude[1]), denscomp_is, rolen-1);
	denscomp_is[rolen-1] = 2*gradientMagnitude[rolen-1];
	scalefloats(denscomp_is, rolen, 0.5);

	multiplyfloats(denscomp_is, denscomp_it, denscomp, rolen);

	/*Deallocate mem*/
	free(kr);
	free(gradientMagnitude);
	free(Gtwist);
	free(denscomp_it);
	free(denscomp_is);
	free(kPhasexy);
	free(kMagxy);

	return;
}


struct Trajectory* generateSpirals(struct VariableDensity *variableDensity, float fieldOfView, float spatialResolution, float readoutDuration, int rewindTrajectory, float samplingInterval, int interleavesDesired, enum SpiralType sptype, float floretAngle, float readoutFieldOfView, float maxGradientAmplitude, float maxSlewRate)
{
  struct Trajectory *trajectory = newTrajectory();
	adjustSpatialResolution(fieldOfView, trajectory->imageDimensions, &spatialResolution);
	float kr = 0;
	float kSpaceMaxRadial = 5/spatialResolution;

	int n;
	int readoutPointsTest;

	if(readoutFieldOfView)
		trajectory->maxReadoutGradientAmplitude = fmin(calculateMaxReadoutGradientAmplitude(readoutFieldOfView, samplingInterval), maxGradientAmplitude);
    else
		trajectory->maxReadoutGradientAmplitude = maxGradientAmplitude;

  int numReadoutPointsDesired = readoutDuration/samplingInterval;
  numReadoutPointsDesired += numReadoutPointsDesired%2;

	int interleavesLow = 1;
	int interleavesHigh = 0;

	int steps = 1;
	if(variableDensity)
		steps = variableDensity->steps;

	   for(n=0; n<steps; n++)
	   {
		   float fieldOfViewFinal;
		   if(variableDensity)
		   {
				kr = variableDensity->step[n].kr;
			getFinalFieldOfView(variableDensity, &fieldOfView, &fieldOfViewFinal, 1);
		   }
		   else
		   {
			   fieldOfViewFinal = fieldOfView;
			   kr = kSpaceMaxRadial;
		   }

		   if(sptype==spFERMAT)
          fieldOfViewFinal = calcFovFermatFloret(fieldOfViewFinal, kr, floretAngle);

			interleavesHigh = fmax(interleavesHigh, 2.0f*M_PI*kr*fieldOfViewFinal);
	   }

	   if(kr<kSpaceMaxRadial && variableDensity)
       {
		   float fieldOfViewAtKr = fieldOfView;
			getFinalFieldOfView(variableDensity, &fieldOfView, &fieldOfViewAtKr, 1);
          if(sptype==spFERMAT)
            fieldOfViewAtKr = calcFovFermatFloret(fieldOfViewAtKr, kSpaceMaxRadial, floretAngle);

		   interleavesHigh = fmax(interleavesHigh, 2.0f*M_PI*kr*fieldOfViewAtKr);
       }

    float* gradientReadoutWaveformsBasis = (float*)calloc(2*numReadoutPointsDesired, sizeof(float));
		float* gradientWaveformsTest = NULL;
		float* kSpaceCoordinatesTest = NULL;

		while((interleavesHigh-interleavesLow)>1)
        {
			int interleavesTest = (interleavesLow+interleavesHigh)/2.0f;

			if(!readoutDuration)
			{
				interleavesTest = interleavesDesired;
        trajectory->numReadoutPoints = 1e5;
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

      if(!generateSpiral(fieldOfView, spatialResolution, variableDensity, interleavesTest, numReadoutPointsDesired, samplingInterval, sptype, floretAngle, trajectory->maxReadoutGradientAmplitude, maxSlewRate, &kSpaceCoordinatesTest, &gradientWaveformsTest, &readoutPointsTest))
            {
        if(readoutPointsTest > numReadoutPointsDesired)
					interleavesLow = interleavesTest;
                else
                {
					interleavesHigh = interleavesTest;
          trajectory->numReadoutPoints = readoutPointsTest;
          memcpy(gradientReadoutWaveformsBasis, gradientWaveformsTest, 2*trajectory->numReadoutPoints*sizeof(float));
          trajectory->numReadouts = interleavesTest;
                }
            }
            else
				interleavesLow = interleavesTest;

			printf(" %d interleaves %d pts\n", interleavesTest, readoutPointsTest);

			if(!readoutDuration)
				interleavesHigh = interleavesLow;
        }

    float* gradientWaveformsBasis[2];
    if(rewindTrajectory)
    {
      traverseKspaceToZero(gradientReadoutWaveformsBasis, &gradientReadoutWaveformsBasis[trajectory->numReadoutPoints], NULL, trajectory->numReadoutPoints, samplingInterval, maxGradientAmplitude, maxSlewRate, &gradientWaveformsBasis[0], &gradientWaveformsBasis[1], NULL, &trajectory->numWaveformPoints);
    }
    else
    {
      gradientWaveformsBasis[0] = gradientReadoutWaveformsBasis;
      gradientWaveformsBasis[1] = &gradientReadoutWaveformsBasis[trajectory->numReadoutPoints];
      trajectory->numWaveformPoints = trajectory->numReadoutPoints;
    }

  allocateTrajectory(trajectory, trajectory->numReadoutPoints, trajectory->numWaveformPoints, 2, 1, trajectory->numReadouts, STORE_ALL);
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
  trajectory->numPreReadoutPoints = 0;
	rotateBasis(gradientWaveformsBasis[0], gradientWaveformsBasis[1], trajectory, 2*M_PI);

  for(n=0; n<trajectory->numReadouts; n++)
	{
    float* gx = &trajectory->gradientWaveforms[2*n*trajectory->numWaveformPoints];
    float* gy = &trajectory->gradientWaveforms[(2*n+1)*trajectory->numWaveformPoints];

    float* kx = &trajectory->kSpaceCoordinates[2*n*trajectory->numReadoutPoints];
    float* ky = &trajectory->kSpaceCoordinates[(2*n+1)*trajectory->numReadoutPoints];
    calcSpiralDcf(gx, gy, kx, ky, trajectory->numReadoutPoints, &trajectory->densityCompensation[n*trajectory->numReadoutPoints]);
	}

    return trajectory;
}


