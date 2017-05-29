/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#include "cones.h"

#include "trajectory.h"
#include "angles.h"
#include "mathops.h"
#include "arrayops.h"
#include "variabledensity.h"
#include "convertendian.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

const float DESIGN_SAMPLING_INTERVAL = 1e-6;

#ifndef GYROMAGNETIC_RATIO
#define GYROMAGNETIC_RATIO 4257
#endif

double calcc1(double sintheta, double fieldOfViewRadial, double interleaves)
{
	return 4*M_PI*M_PI*sintheta*sintheta*fieldOfViewRadial*fieldOfViewRadial/(interleaves*interleaves);
}

double calcc2(double fieldOfViewCircumferential, double fieldOfViewRadial)
{
	return fieldOfViewRadial*fieldOfViewRadial/(fieldOfViewCircumferential*fieldOfViewCircumferential);
}

float calculateInterconeCompensation(enum InterConeCompensation interconeCompensation, float kr, float kSpaceExtentRadial)
{
	float compensation;

	switch(interconeCompensation)
	{
		case Compensation1:
			compensation = 1 + sqrt(fmax(0, 1-4*kr*kr/(kSpaceExtentRadial*kSpaceExtentRadial)));
			break;
		default:
			compensation = 1.0f;
			break;
	}

	return compensation;
}

float calculateGtwist(float kr, float kSpaceExtentRadial, float theta, float fieldOfViewRadial, float fieldOfViewCircumferential, float interleaves, enum InterConeCompensation interconeCompensation)
{
	float gtwist;
	float c1 = calcc1(sin(theta), fieldOfViewRadial, interleaves);
	float c2 = calcc2(fieldOfViewCircumferential, fieldOfViewRadial);
	float compensation = 1.0f;

//	if(ic==1)
//		icomp = calcConesComp(ic, kSpaceRadius, kmax);

//	if(ic>1)
	compensation = calculateInterconeCompensation(interconeCompensation, kr, kSpaceExtentRadial);
//		c2 = 0;
	c2 = 0;

	gtwist = c1*kr*kr/(compensation*compensation) - c2;
	gtwist = sqrt((gtwist>0)?gtwist:0);

	return gtwist;
}

int generateCone(float fieldOfViewRadial, float fieldOfViewCircumferential, const struct VariableDensity *variableDensity, float kSpaceExtentRadial, float interleaves, enum InterConeCompensation compensation, float theta, int maxPoints, float samplingInterval, int rotatable, float scaleXY, float scaleZ, float maxGradientAmplitude, float maxSlewRate, int *readoutPoints, float **kSpaceCoordinates, float **gradientWaveforms)
{
	int dimensions = 3;
	int oversamplingRatio = roundToInteger(samplingInterval/DESIGN_SAMPLING_INTERVAL);
	float gammaDeltaTime = GYROMAGNETIC_RATIO*DESIGN_SAMPLING_INTERVAL;

	float maxGradientAmplitudeXY = scaleXY*maxGradientAmplitude;
	float maxGradientAmplitudeZ = scaleZ*maxGradientAmplitude;
	maxSlewRate *= .996;
	float maxSlewRateXY = scaleXY*maxSlewRate;
	float maxSlewRateZ = scaleZ*maxSlewRate;
	float maxGradientDeltaXY = maxSlewRateXY*DESIGN_SAMPLING_INTERVAL;
	float maxGradientDeltaZ = maxSlewRateZ*DESIGN_SAMPLING_INTERVAL;
	float maxGradientDelta = maxSlewRate*DESIGN_SAMPLING_INTERVAL;

	int maxDesignPoints = oversamplingRatio*maxPoints;
	float *kDesign = (float*)calloc(dimensions*maxDesignPoints, sizeof(float));
	int *stepSign = (int*)malloc(maxDesignPoints*sizeof(int));
	//float *twist = (float*)malloc(maxDesignPoints*sizeof(float));
	//twist[0] = 0;
	float twist;
	float scale = 1;

	int n;
	for(n=0; n<maxDesignPoints; n++)
		stepSign[n] = 1;

	/* set transverse and longitudinal kmaxs */
	int d;
//	for(d=0; d<2; d++)
//		kmaxs[d] = 5/res[d];

	/*theta = PI/2.0f - theta;*/
	float costheta = cos(theta);
	float sintheta = sin(theta);

	/* calculate kmax at theta */
//	kmax = intlineellipse(kmaxs[1], kmaxs[0], theta);

//	fovc = vd->fov[0];
//	fovr = 1/intlineellipse(1.0f/vd->fov[1], 1.0f/vd->fov[0], theta);
	/*gtwiststart = calcgtwist(*kstart, kmax, theta, fovr, fovc, interleaves, ic);*/

	kDesign[0] = gammaDeltaTime*maxGradientDelta*sintheta;
	kDesign[2] = gammaDeltaTime*maxGradientDelta*costheta;

	kDesign[3] = 2*gammaDeltaTime*maxGradientDelta*sintheta;
	kDesign[5] = 2*gammaDeltaTime*maxGradientDelta*costheta;

	float kr = 2*gammaDeltaTime*maxGradientDelta;

	n = 1;
	int iteration = 1;
//	*ntwist = -1;
	int stepBack;
	float gradientMagnitude=0;
	int valid = 1;
	float fieldOfViewInitial[2] = {fieldOfViewRadial, fieldOfViewCircumferential};
	float fieldOfViewCurrent[2];
	float kDesignNext[3];
	float gradient[3];
	float deltaK[3];
	float kSpaceDirection[3];
	float kSpaceDirectionMagnitudeXY;
	float projectedGradientMagnitude;
	float projectedGradientMagnitudeXY;
	float term;
	float gradientMagnitudeRange[2];
	float gradientMagnitudeRangeXY[2];
	float gradientMagnitudeRangeZ[2];
	float gradientMagnitudeXY;
	float gradientMagnitudeZ;
	float s;
	float phi;
	float deltaPhiDeltaKr;
	float deltaKr;
	float deltaPhi;
	float phiNext;
	float krNext;

	while((kr <= kSpaceExtentRadial+oversamplingRatio*gradientMagnitude*gammaDeltaTime) && valid)
	{
		stepBack = 0;

//		for(d=0; d<dimensions; d++)
//			km[d] = 1.5*kDesign[n*dimensions+d] - 0.5*kDesign[(n-1)*dimensions+d];
		/*kmr = norm2(km, dimensions);*/
		/*knorm = kmr/kmax;*/
		if(variableDensity)
			getFieldOfView(variableDensity, kr, fieldOfViewInitial, fieldOfViewCurrent, 2);
		else
			memcpy(fieldOfViewCurrent, fieldOfViewInitial, 2*sizeof(float));

//		fovc = fov[0];
//		fovr = 1/intlineellipse(1/fov[1], 1/fov[0], theta);

		twist = calculateGtwist(kr, kSpaceExtentRadial, theta, fieldOfViewRadial, fieldOfViewCircumferential, interleaves, compensation);

		/*if(ic==2 && gtwist)
		{
			projectedGradient = gtwist*gtwist-calcc2(vd->fov[0], fovr);
			gtwist = sqrt((projectedGradient>0)?projectedGradient:0);
		}*/

//		twist[n] = twistCurrent;

		phi = atan2(kDesign[n*dimensions+1], kDesign[n*dimensions]);
		deltaPhiDeltaKr = twist/(kr*sintheta);
		/*deltaKspaceRadius = kmr-kSpaceRadius;*/
		deltaKr = norm2(&kDesign[n*dimensions], dimensions) - norm2(&kDesign[(n-1)*dimensions], dimensions);
		deltaPhi = deltaPhiDeltaKr*deltaKr;
		phiNext = phi+deltaPhi;

		kDesignNext[0] = cos(phiNext)*sintheta;
		kDesignNext[1] = sin(phiNext)*sintheta;
		kDesignNext[2] = costheta;
		krNext = kr+deltaKr;
		scalefloats(kDesignNext, dimensions, krNext);
		/*scalefloats(knext, dimensions, kSpaceRadius+deltaKspaceRadius);*/
		subtractArray(kDesignNext, &kDesign[n*dimensions], deltaK, dimensions);
//		copyfloats(deltaK, kSpaceDirection, dimensions);
		memcpy(kSpaceDirection, deltaK, dimensions*sizeof(float));
		scalefloats(kSpaceDirection, dimensions, 1.0f/norm2(deltaK,dimensions));

		for(d=0; d<dimensions; d++)
			gradient[d] = (kDesign[n*dimensions+d]-kDesign[(n-1)*dimensions+d])/gammaDeltaTime;

		if(rotatable)
		{
			/*projectedGradient = u[0]*g[0]/scaleXY + u[1]*g[1]/scaleXY + u[2]*g[2]/scaleZ;
			term = (maxGradientDelta*maxGradientDelta - ((g[0]*g[0] + g[1]*g[1])/(scaleXY*scaleXY) + g[2]*g[2]/(scaleZ*scaleZ))) + projectedGradient*projectedGradient;*/
				projectedGradientMagnitude = dot(kSpaceDirection, gradient, dimensions);
				term = (maxGradientDelta*maxGradientDelta - (gradient[0]*gradient[0] + gradient[1]*gradient[1] + gradient[2]*gradient[2])) + projectedGradientMagnitude*projectedGradientMagnitude;
		}
		else
		{
			kSpaceDirectionMagnitudeXY = sqrt(kSpaceDirection[0]*kSpaceDirection[0]+kSpaceDirection[1]*kSpaceDirection[1]);
			projectedGradientMagnitudeXY = dot(kSpaceDirection, gradient, 2);
			//float projectedGradientMagnitudeZ = kSpaceDirection[2]*gradient[2];

			/*termscaleZ = u[2]*u[2]*(maxGradientDelta*maxGradientDelta - g[2]*g[2]/(scaleZ*scaleZ)) + projectedGradientz*projectedGradientz/(scaleZ*scaleZ);*/
			/*if(termscaleXY>=0 && termscaleZ>=0)*/
				term = kSpaceDirectionMagnitudeXY*kSpaceDirectionMagnitudeXY*(maxGradientDeltaXY*maxGradientDeltaXY - (gradient[0]*gradient[0] + gradient[1]*gradient[1])) + projectedGradientMagnitudeXY*projectedGradientMagnitudeXY;
			/*else
				term = -1;*/
		}

		if(term>=0)
		{
			for(d=0; d<2; d++)
			{
				if(d==0)
					s = -1.0f;
				else
					s = 1.0f;

				if(rotatable)
				{
					gradientMagnitudeRange[d] = (projectedGradientMagnitude + s*sqrt(term));
				}
				else
				{
					gradientMagnitudeRangeXY[d] = (projectedGradientMagnitudeXY + s*sqrt(term))/(kSpaceDirectionMagnitudeXY*kSpaceDirectionMagnitudeXY);
					/*gradientMagnitudeRangeZ[d] = (g[2] + gs*maxGradientDelta)/fabsf(u[2]);*/
					gradientMagnitudeRangeZ[d] = (gradient[2] + s*maxGradientDeltaZ)/fabs(kSpaceDirection[2]);
				}
			}

			if(rotatable)
			{
				gradientMagnitudeRange[1] = fmin(gradientMagnitudeRange[1], maxGradientAmplitude);
				gradientMagnitudeRange[1] = fmax(gradientMagnitudeRange[1], -maxGradientAmplitude);
			}
			else
			{
				gradientMagnitudeRange[0] = fmax(gradientMagnitudeRangeXY[0],gradientMagnitudeRangeZ[0]);
				gradientMagnitudeRange[1] = fmin(gradientMagnitudeRangeXY[1],gradientMagnitudeRangeZ[1]);
			}

			if(gradientMagnitudeRange[0]>gradientMagnitudeRange[1])
				stepBack = 1;
		}
		else
			stepBack = 1;

		if(rotatable && gradientMagnitudeRange[0]>maxGradientAmplitude)
			stepBack = 1;

		if(!stepBack)
		{
			if(stepSign[n]==1)
				gradientMagnitude = gradientMagnitudeRange[1];
			else
				gradientMagnitude = gradientMagnitudeRange[0];

			if(!rotatable)
			{
				/*gradientMagnitudeXY = gm*kSpaceDirectionMagnitudeXY/scaleXY;
				gradientMagnitudeZ = gm*u[2]/scaleZ;*/
				gradientMagnitudeXY = gradientMagnitude*kSpaceDirectionMagnitudeXY;
				gradientMagnitudeZ = gradientMagnitude*kSpaceDirection[2];


				scale = 1.0f;
				if(gradientMagnitudeXY>maxGradientAmplitudeXY)
					scale = maxGradientAmplitudeXY/gradientMagnitudeXY;

				if(gradientMagnitudeZ>maxGradientAmplitudeZ)
					scale = fmin(scale, maxGradientAmplitudeZ/gradientMagnitudeZ);

				if(gradientMagnitude>maxGradientAmplitude)
					scale = fmin(scale, maxGradientAmplitude/gradientMagnitude);

				gradientMagnitude *= scale;
				
			}

			for(d=0; d<dimensions; d++)
			{
				gradient[d] = gradientMagnitude*kSpaceDirection[d];
				kDesign[(n+1)*dimensions+d] = kDesign[n*dimensions+d] + gradient[d]*gammaDeltaTime;
			}

//			if(gtwist && *ntwist==-1)
//			{
//				*ntwist = n/oversamplingRatio;
//				*kstart = norm2(&k[n*dimensions], dimensions);
//			}

			n++;
		}
		else
		{
			while((n>4) && (stepSign[n-1]==-1))
				n--;
			stepSign[n-1] = -1;
			n -= 2;
		}
		kr = norm2(&kDesign[n*dimensions], dimensions);

		valid = n<(maxDesignPoints-1);
		iteration++;
	}

	if(valid)
	{
		*readoutPoints = (n-1)/oversamplingRatio;

		*kSpaceCoordinates = (float*)malloc(*readoutPoints*dimensions*sizeof(float));
		*gradientWaveforms = (float*)malloc(*readoutPoints*dimensions*sizeof(float));
//		*Gtwistout = (float*)malloc(*readoutPoints*sizeof(float));

		/*for(d=0; d<dimensions; d++)
		{
			(*kSpaceCoordinates)[*Nk*d] = k[d];
			(*gradientWaveforms)[*Nk*d] = k[d]/gammaDeltaTime;

			if(d<3)
			{
				(*kSpaceCoordinates)[*Nk*d] /= scaleXY;
				(*gradientWaveforms)[*Nk*d] /= scaleXY;
			}
			else
			{
				(*kSpaceCoordinates)[*Nk*d] /= scaleZ;
				(*gradientWaveforms)[*Nk*d] /= scaleZ;
			}
		}

		(*Gtwistout)[0] = 0;*/
		for(n=0; n<*readoutPoints; n++)
		{
//			(*Gtwistout)[n] = Gtwist[oversamplingRatio*n];
			for(d=0; d<dimensions; d++)
			{
				(*kSpaceCoordinates)[*readoutPoints*d+n] = kDesign[oversamplingRatio*n*dimensions+d];
				if(n)
					(*gradientWaveforms)[*readoutPoints*d+n] = (kDesign[oversamplingRatio*n*dimensions+d]-kDesign[oversamplingRatio*(n-1)*dimensions+d])/gammaDeltaTime/oversamplingRatio;
				else
					(*gradientWaveforms)[*readoutPoints*d+n] = kDesign[oversamplingRatio*n*dimensions+d]/gammaDeltaTime/oversamplingRatio;

			}
		}
	}

	free(kDesign);
	free(stepSign);
//	free(Gtwist);

	return valid;
}

void interpolateCone(float fieldOfViewXY, float fieldOfViewZ, float spatialResolutionXY, float spatialResolutionZ, const float *basisGradientWaveforms, int basisWaveformPoints, int readoutPoints, int waveformPoints, float theta, float scaleXY, float scaleZ, struct VariableDensity *variableDensity, int interconeCompensation, float dtheta_ic, float samplingInterval, int* interleaves, float* denscomp, float *interpolatedGradientWaveforms)
{
	float *kSpaceCoordinates;
	float* gradient[3];
	float* kSpaceCoordinate[3];
	int d;
	float *kr, *gradientMagnitude;
	float* Gtwist;
	float* denscomp_ic, *denscomp_it, *denscomp_is;
	float *interleavesAtKr;
	
	float kSpaceExtentRadial;
	float fieldOfViewRadial;
	float fieldOfViewInitial[2] = {fieldOfViewZ, fieldOfViewXY};
	float fieldOfViewAtKr[2];
	float spatialResolution[2] = {spatialResolutionZ, spatialResolutionXY};
	float *kPhaseXY, *kMagXY;
	enum AngleShape shapeFieldOfView = InverseEllipticalShape;
	float kSpaceExtentHalf[2];
	int n;

	kSpaceCoordinates = (float*)malloc(3*readoutPoints*sizeof(float));
	for(d=0; d<3; d++)
	{
		gradient[d] = &interpolatedGradientWaveforms[d*waveformPoints];
		memcpy(gradient[d], &basisGradientWaveforms[d*waveformPoints], basisWaveformPoints*sizeof(float));
		kSpaceCoordinate[d] = &kSpaceCoordinates[d*readoutPoints];
	}

	scalefloats(gradient[0], basisWaveformPoints, scaleXY);
	scalefloats(gradient[1], basisWaveformPoints, scaleXY);
	scalefloats(gradient[2], basisWaveformPoints, scaleZ);

	for(d=0; d<3; d++)
		gradientToKspace(gradient[d], kSpaceCoordinate[d], samplingInterval, readoutPoints);

	kr = (float*)malloc(readoutPoints*sizeof(float));
	gradientMagnitude = (float*)malloc(readoutPoints*sizeof(float));
	for(n=0; n<readoutPoints; n++)
	{
		kr[n] = sqrt(kSpaceCoordinate[0][n]*kSpaceCoordinate[0][n] + kSpaceCoordinate[1][n]*kSpaceCoordinate[1][n] + kSpaceCoordinate[2][n]*kSpaceCoordinate[2][n]);
		gradientMagnitude[n] = sqrt(gradient[0][n]*gradient[0][n] + gradient[1][n]*gradient[1][n] + gradient[2][n]*gradient[2][n]);
	}

	Gtwist = (float*)malloc(readoutPoints*sizeof(float));
	kPhaseXY = (float*)malloc(readoutPoints*sizeof(float));
	kMagXY = (float*)malloc(readoutPoints*sizeof(float));

	calculatePhase(kSpaceCoordinate[0], kSpaceCoordinate[1], kPhaseXY, readoutPoints, 0, 1.0);
	calculateMagnitude(kSpaceCoordinate[0], kSpaceCoordinate[1], kMagXY, readoutPoints);
	unwrapPhase(kPhaseXY, Gtwist, readoutPoints, readoutPoints, M_PI);
	for(n=0; n<readoutPoints-1; n++)
        Gtwist[n] = fmax((Gtwist[n+1]-Gtwist[n])/(kr[n+1]-kr[n])*kMagXY[n], 0);

	interleavesAtKr = (float*)malloc(readoutPoints*sizeof(float));

	/* storing kmaxs */
	for(n=0; n<2; n++)
		kSpaceExtentHalf[n] = 5/spatialResolution[n];

    /*kmax = intlineellipse(5/res[0], 5/res[1], theta);*/
	kSpaceExtentRadial = getExtent(EllipticalShape, theta, kSpaceExtentHalf);

//    if(interconeCompensation>1)
//    {
//		/*FOVrad = 1/intlineellipse(1.0f/vd->fov[0], 1.0f/vd->fov[1], theta);
//		FOVphi = 1/intlineellipse(1.0f/vd->fov[0], 1.0f/vd->fov[1], theta+PI/2.0f);
//	    FOVrad = intlineellipse(vd->fov[0], vd->fov[1], theta);
//	    FOVphi = intlineellipse(vd->fov[0], vd->fov[1], theta+PI/2.0f);
//
//		 might want to use kmax from trajectory rather than theoretical
//		resr = 5/intlineellipse(5/res[0], 5/res[1], theta);*/
//
//		fieldOfViewRadial = getExtent(shapeFieldOfView, theta, fieldOfViewInitial);
//		//FOVtheta = getExtent(shapeFieldOfView, theta+PI/2, vd->fov);
//
//	    /*calcAngleComp(vd, FOVrad, resr);*/
//		/*calcAngleCompa(vd, FOVrad, FOVtheta, kmax);*/
//		calculateInterconeDensity(variableDensity, theta, fieldOfViewXY, fieldOfViewZ, 5/kSpaceExtentRadial, interconeCompensation==3);
//    }
//
    for(n=0; n<readoutPoints-1; n++)
	{
		if(variableDensity)
			getFieldOfView(variableDensity, kr[n], fieldOfViewInitial, fieldOfViewAtKr, 2);
		else
			memcpy(fieldOfViewAtKr, fieldOfViewInitial, 2*sizeof(float));
		/*FOVrad = 1/intlineellipse(1.0f/fov[0], 1.0f/fov[1], theta);*/
		/*FOVrad = intlineellipse(fov[0], fov[1], theta);*/
		fieldOfViewRadial = getExtent(shapeFieldOfView, theta, fieldOfViewAtKr);
		interleavesAtKr[n] = 2*M_PI*kMagXY[n]*fieldOfViewRadial/sqrt(Gtwist[n]*Gtwist[n]+(fieldOfViewRadial/fieldOfViewXY)*(fieldOfViewRadial/fieldOfViewXY));
		switch(interconeCompensation)
		{
			case 1:
				/*nints[n] /= (1+sqrt(max(0,1-4*(kr[n]/kmax)*(kr[n]/kmax))));*/
				interleavesAtKr[n] /= calculateInterconeCompensation(interconeCompensation, kr[n], kSpaceExtentRadial);
				break;
			case 2:
			/*nints[n] = 2*PI*kMagXY[n]*FOVrad/sqrt(Gtwist[n]*Gtwist[n] + (FOVrad/vd->fov[0])*(FOVrad/vd->fov[0]));*/
				break;
		}
	}
		
	*interleaves = ceil(interleavesAtKr[0]);
    for(n=1; n<readoutPoints-1; n++)
	   *interleaves = fmax(*interleaves, ceil(interleavesAtKr[n]));
	/*(*nint)++;*/
	
	//*nint += *nint%2;

	Gtwist[readoutPoints-1] = Gtwist[readoutPoints-2];

	/* Density compensation due to inter-cone spacing
	denscomp_ic = dtheta_ic*kr; */
	denscomp_ic = (float*)malloc(readoutPoints*sizeof(float));
	memcpy(denscomp_ic, kr, readoutPoints*sizeof(float));
	scalefloats(denscomp_ic, readoutPoints, dtheta_ic);

	/* Density compensation due to inter-trajectory spacing (ignoring NINT)
	denscomp_it = abs(kcxy)./sqrt(1+Gtwist.^2);*/
	denscomp_it = (float*)malloc(readoutPoints*sizeof(float));
	for(n=0; n<readoutPoints; n++)
		denscomp_it[n] = kMagXY[n]/sqrt(1+Gtwist[n]*Gtwist[n]);

	/* Density compensation due to inter-sample spacing
	denscomp_is = (gr+[gr(2:end); gr(end)])/2;*/
	denscomp_is = (float*)malloc(readoutPoints*sizeof(float));
	addfloats(gradientMagnitude, &(gradientMagnitude[1]), denscomp_is, readoutPoints-1);
	denscomp_is[readoutPoints-1] = 2*gradientMagnitude[readoutPoints-1];
	scalefloats(denscomp_is, readoutPoints, 0.5);

	if(kMagXY[readoutPoints-1]/kr[readoutPoints-1] < 1e-8)
	{
		/* Special Case of the null cone
		denscomp_ic = (dtheta_ic^2*(kr).^2)/8;*/
		multiplyfloats(kr, kr, denscomp_ic, readoutPoints);
		scalefloats(denscomp_ic, readoutPoints, dtheta_ic*dtheta_ic/8);

	}

	multiplyfloats(denscomp_ic, denscomp_it, denscomp, readoutPoints);
    multiplyfloats(denscomp, denscomp_is, denscomp, readoutPoints);

	free(kSpaceCoordinates);
	free(kr);
	free(gradientMagnitude);
	free(Gtwist);
	free(interleavesAtKr);
	free(denscomp_ic);
	free(denscomp_it);
	free(denscomp_is);
	free(kPhaseXY);
	free(kMagXY);

	return;
}

void initializeConesInterpolation(struct ConesInterpolation *interpolation, int readouts)
{
	interpolation->readouts = readouts;
	interpolation->theta = (float*)malloc(readouts*sizeof(float));
	interpolation->thetaIndex = (int*)malloc(readouts*sizeof(int));
	interpolation->scaleXY = (float*)malloc(readouts*sizeof(float));
	interpolation->scaleZ = (float*)malloc(readouts*sizeof(float));
	interpolation->phi = (float*)malloc(readouts*sizeof(float));
	interpolation->cone = (int*)malloc(readouts*sizeof(int));
	interpolation->readout = (int*)malloc(readouts*sizeof(int));
	interpolation->interleavesOnCone = (int*)malloc(readouts*sizeof(int));
	interpolation->interleafOnCone = (int*)malloc(readouts*sizeof(int));
}

void makeConesInterpolation(struct Cones *cones)
{
	struct ConesInterpolation singleInterleafInterpolation;
	struct Trajectory *trajectory = &cones->trajectory;
	float thetaMax = M_PI_2;
	float coneCoverageStart;
	float coneCoverageEnd;
	float deltaConeAngle;
	const float *basisGradientWaveforms;
	float *interpolatedGradientWaveforms = (float*)malloc(3*trajectory->waveformPoints*sizeof(float));
	float phiOffset;

	int c;
	int r=0;
	int i;
	initializeConesInterpolation(&singleInterleafInterpolation, cones->numCones);
	trajectory->readouts = 0;
	cones->coneAngleDensityCompensation = (float*)malloc(cones->numCones*trajectory->readoutPoints*sizeof(float));
	for (c=0; c<cones->numCones; c++)
	{
		singleInterleafInterpolation.theta[c] = cones->coneAngles[c];
//		singleInterleafSchedule.coneIndex[c] = fmin(trajectory->bases, fmax(1,ceil(fabs(singleInterleafSchedule.theta[c])/(M_PI_2)*trajectory->bases)))-1;
		singleInterleafInterpolation.cone[c] = (1-fabs(1-singleInterleafInterpolation.theta[c]*M_2_PI))*trajectory->bases;

		coneCoverageStart = singleInterleafInterpolation.cone[c]/(1.0*trajectory->bases)*thetaMax;
		singleInterleafInterpolation.scaleXY[c] = cos(singleInterleafInterpolation.theta[c])/cos(coneCoverageStart);

		coneCoverageEnd = (singleInterleafInterpolation.cone[c]+1)/(1.0*trajectory->bases)*thetaMax;
		singleInterleafInterpolation.scaleZ[c] = sin(singleInterleafInterpolation.theta[c])/sin(coneCoverageEnd);

	   if(c)
		   deltaConeAngle = singleInterleafInterpolation.theta[c]-singleInterleafInterpolation.theta[c-1];
	   else
		   deltaConeAngle = 2*singleInterleafInterpolation.theta[c];

		basisGradientWaveforms = &cones->basisGradientWaveforms[3*trajectory->waveformPoints*singleInterleafInterpolation.cone[c]];
		interpolateCone(trajectory->fieldOfView[0], trajectory->fieldOfView[2], trajectory->spatialResolution[0], trajectory->spatialResolution[2], basisGradientWaveforms, cones->basisWaveformPoints[singleInterleafInterpolation.cone[c]], cones->basisReadoutPoints[singleInterleafInterpolation.cone[c]], trajectory->waveformPoints, singleInterleafInterpolation.theta[c], singleInterleafInterpolation.scaleXY[c], singleInterleafInterpolation.scaleZ[c], trajectory->variableDensity, cones->interconeCompensation, deltaConeAngle, trajectory->samplingInterval, &singleInterleafInterpolation.interleavesOnCone[c], &cones->coneAngleDensityCompensation[c*trajectory->readoutPoints], interpolatedGradientWaveforms);
		trajectory->readouts += singleInterleafInterpolation.interleavesOnCone[c];

		printf("theta %f\tcone %d\tnintl %d\n", singleInterleafInterpolation.theta[c], singleInterleafInterpolation.cone[c], singleInterleafInterpolation.interleavesOnCone[c]);
	}
	printf("Max readout length:\t%d\n", trajectory->readoutPoints);
	printf("readouts %d\n", trajectory->readouts);

	initializeConesInterpolation(&cones->interpolation, trajectory->readouts);
	for(c=0; c<cones->numCones; c++)
	{
		if(c)
			phiOffset = M_PI/singleInterleafInterpolation.interleavesOnCone[c-1];
		else
			phiOffset = 0;

		for(i=0; i<singleInterleafInterpolation.interleavesOnCone[c]; i++)
		{
			cones->interpolation.readout[r] = r;
			cones->interpolation.cone[r] = singleInterleafInterpolation.cone[c];
			cones->interpolation.scaleXY[r] = singleInterleafInterpolation.scaleXY[c];
			cones->interpolation.scaleZ[r] = singleInterleafInterpolation.scaleZ[c];
			cones->interpolation.theta[r] = singleInterleafInterpolation.theta[c];
			cones->interpolation.thetaIndex[r] = c;
			cones->interpolation.interleavesOnCone[r] = singleInterleafInterpolation.interleavesOnCone[c];
			cones->interpolation.phi[r] = i*2*M_PI/singleInterleafInterpolation.interleavesOnCone[c] + phiOffset;
		}
	}
}


int generateConesBasis(struct Cones *cones)
{
	int status = 0;
	enum AngleShape elevationAngleFieldOfViewShape = InverseEllipticalShape;

	float kSpaceExtent[2];
	float initialElevationAngleDelta;
	struct Trajectory *trajectory = &cones->trajectory;

	struct VariableDensity *variableDensity = cones->trajectory.variableDensity;

	float fieldOfView[2] = {trajectory->fieldOfView[2], trajectory->fieldOfView[0]};

	float finalFieldOfView[2];

	int d;
	int b;
	float fieldOfViewRadial;
	float fieldOfViewCircumferential;
	float spatialResolutionRadial;
	float kSpaceExtentRadial;
	float elevationAngleParametric;
	float scaleXY;
	float scaleZ;
	float interleavesLow;
	float interleavesHigh;
	float interleaves;
	float *basisReadoutGradientWaveforms;
	int currentReadoutPoints;
	float *currentKspaceCoordinates = NULL;
	float *currentGradientWaveforms = NULL;
	float **basisGradientRewoundX;
	float **basisGradientRewoundY;
	float **basisGradientRewoundZ;

	printf("Number of readout points:\t%d\n", trajectory->readoutPoints);
	printf("Number of basis waveforms:\t%d\n", trajectory->bases);

	srand(0);

	for(d=0; d<2; d++)
		kSpaceExtent[d] = 10/trajectory->spatialResolution[2-d];
	initialElevationAngleDelta = 1/(kSpaceExtent[0]*fieldOfView[0]);

	if(variableDensity)
		getFinalFieldOfView(variableDensity, fieldOfView, finalFieldOfView, 2);
	else
		memcpy(finalFieldOfView, fieldOfView, 2*sizeof(float));

	calculateAngles(initialElevationAngleDelta/2, M_PI-initialElevationAngleDelta/2, elevationAngleFieldOfViewShape, finalFieldOfView, EllipticalShape, kSpaceExtent, &cones->coneAngles, NULL, NULL, &cones->numCones);

	printf("Number of cones:\t%d\n", cones->numCones);

	cones->designConeAngles = (float*)malloc(trajectory->bases*sizeof(float));

	basisReadoutGradientWaveforms = (float*)malloc(3*trajectory->bases*trajectory->readoutPoints*sizeof(float));
	
	cones->basisKspaceCoordinates = (float*)malloc(3*trajectory->bases*trajectory->readoutPoints*sizeof(float));
	basisGradientRewoundX = (float**)malloc(trajectory->bases*sizeof(float*));
	basisGradientRewoundY = (float**)malloc(trajectory->bases*sizeof(float*));
	basisGradientRewoundZ = (float**)malloc(trajectory->bases*sizeof(float*));
	for(b=0; b<trajectory->bases; b++)
	{
		printf("Waveform %d\n", b+1);
		float fromElevationAngle = M_PI_2*b/trajectory->bases;
		float toElevationAngle = M_PI_2*(b+1.0f)/trajectory->bases;

		cones->designConeAngles[b] = atan2(kSpaceExtent[0]*sin(toElevationAngle), kSpaceExtent[1]*cos(fromElevationAngle));

		fieldOfViewRadial = getExtent(InverseEllipticalShape, cones->designConeAngles[b], fieldOfView);
		fieldOfViewCircumferential = getExtent(InverseEllipticalShape, cones->designConeAngles[b]+M_PI_2, fieldOfView);
		spatialResolutionRadial = 5/getExtent(EllipticalShape, cones->designConeAngles[b], kSpaceExtent);
		kSpaceExtentRadial = 5/spatialResolutionRadial;
		elevationAngleParametric = atan2((kSpaceExtentRadial*sin(cones->designConeAngles[b])*kSpaceExtent[1]),(kSpaceExtentRadial*cos(cones->designConeAngles[b])*kSpaceExtent[0]));

		scaleXY = cos(elevationAngleParametric)/cos(fromElevationAngle);
		scaleZ = sin(elevationAngleParametric)/sin(toElevationAngle);

		interleavesLow = .01;
		interleavesHigh = 2*M_PI*kSpaceExtentRadial*fieldOfViewRadial*cos(cones->designConeAngles[b]);

		cones->basisReadoutPoints[b] = -1;
		while(interleavesHigh-interleavesLow>.03 && cones->basisReadoutPoints[b]!=trajectory->readoutPoints)
		{
			interleaves = .5*(interleavesLow+interleavesHigh);

			if(currentKspaceCoordinates)
			{
				free(currentKspaceCoordinates);
				currentKspaceCoordinates = NULL;
			}

			if(currentGradientWaveforms)
			{
				free(currentGradientWaveforms);
				currentGradientWaveforms = NULL;
			}

			if(generateCone(fieldOfViewRadial, fieldOfViewCircumferential, variableDensity, kSpaceExtentRadial, interleaves, cones->interconeCompensation, cones->designConeAngles[b], trajectory->readoutPoints, trajectory->samplingInterval, cones->rotatable, scaleXY, scaleZ, trajectory->maxReadoutGradientAmplitude, trajectory->maxSlewRate, &currentReadoutPoints, &currentKspaceCoordinates, &currentGradientWaveforms) && cones->basisReadoutPoints[b]<=trajectory->readoutPoints)
			{
				interleavesHigh = interleaves;
				cones->basisReadoutPoints[b] = currentReadoutPoints;
				for (d=0; d<3; d++)
				{
					memcpy(&basisReadoutGradientWaveforms[(b*3+d)*trajectory->readoutPoints], &currentGradientWaveforms[d*currentReadoutPoints], currentReadoutPoints*sizeof(float));
					memcpy(&cones->basisKspaceCoordinates[(b*3+d)*trajectory->readoutPoints], &currentKspaceCoordinates[d*currentReadoutPoints], currentReadoutPoints*sizeof(float));
				}
			}
			else
			{
				interleavesLow = interleaves;
			}
			printf("Basis %d interleaves %f readout points %d\n", b+1, interleaves, cones->basisReadoutPoints[b]);
		}
	}

	trajectory->waveformPoints = 0;
	for(b=0; b<trajectory->bases; b++)
	{
		/*if(gradientRewoundX)
		{
			free(gradientRewoundX);
			gradientRewoundX = NULL;
		}
		if(gradientRewoundY)
		{
			free(gradientRewoundY);
			gradientRewoundY = NULL;
		}
		if(gradientRewoundZ)
		{
			free(gradientRewoundZ);
			gradientRewoundZ = NULL;
		}*/
traverseKspaceToZero(&basisReadoutGradientWaveforms[b*trajectory->readoutPoints*3], &basisReadoutGradientWaveforms[(b*3+1)*trajectory->readoutPoints], &basisReadoutGradientWaveforms[(b*3+2)*trajectory->readoutPoints], cones->basisReadoutPoints[b], trajectory->samplingInterval, trajectory->maxGradientAmplitude, trajectory->maxSlewRate, &basisGradientRewoundX[b], &basisGradientRewoundY[b], &basisGradientRewoundZ[b], &cones->basisWaveformPoints[b]);
		trajectory->waveformPoints = fmax(trajectory->waveformPoints, cones->basisWaveformPoints[b]);
	}

	trajectory->waveformPoints += trajectory->waveformPoints%2;
	cones->basisGradientWaveforms = (float*)malloc(3*trajectory->bases*trajectory->waveformPoints*sizeof(float));
	for(b=0; b<trajectory->bases; b++)
	{
		memcpy(&cones->basisGradientWaveforms[b*3*trajectory->waveformPoints], basisGradientRewoundX[b], cones->basisWaveformPoints[b]*sizeof(float));
		memcpy(&cones->basisGradientWaveforms[(b*3+1)*trajectory->waveformPoints], basisGradientRewoundY[b], cones->basisWaveformPoints[b]*sizeof(float));
		memcpy(&cones->basisGradientWaveforms[(b*3+2)*trajectory->waveformPoints], basisGradientRewoundZ[b], cones->basisWaveformPoints[b]*sizeof(float));
	}

	return status;
}

float* conesBasisGradientWaveform(const struct Cones* cones, int basis, int axis)
{
	return &cones->basisGradientWaveforms[(basis*3+axis)*cones->trajectory.waveformPoints];
}

void generateReadoutWaveforms(int index, const struct Cones* cones, float* gradientX, float* gradientY, float* gradientZ)
{
	const struct ConesInterpolation *interpolation = &cones->interpolation;
	const struct Trajectory *trajectory = &cones->trajectory;
	int cone = interpolation->cone[index];

	memcpy(gradientX, &cones->basisGradientWaveforms[cone*3*trajectory->waveformPoints], trajectory->waveformPoints*sizeof(float));
	memcpy(gradientY, &cones->basisGradientWaveforms[(cone*3+1)*trajectory->waveformPoints], trajectory->waveformPoints*sizeof(float));
	memcpy(gradientZ, &cones->basisGradientWaveforms[(cone*3+2)*trajectory->waveformPoints], trajectory->waveformPoints*sizeof(float));

	scalecomplex(gradientX, gradientY, interpolation->scaleXY[index]*cos(interpolation->phi[index]), interpolation->scaleXY[index]*sin(interpolation->phi[index]), trajectory->waveformPoints);
	scalefloats(gradientZ, trajectory->waveformPoints, interpolation->scaleZ[index]);
}

struct Cones *generateCones(float fieldOfViewXY, float fieldOfViewZ, const struct VariableDensity *variableDensity, float xySpatialResolution, float zSpatialResolution, int bases, int rotatable, enum InterConeCompensation interConeCompensation, float readoutDuration, float samplingInterval, float filterFieldOfView, float maxGradientAmplitude, float maxSlewRate, enum WaveformStorageType storage)
{
	struct Cones* cones = allocateCones(bases);
	struct Trajectory *trajectory = &cones->trajectory;
	initializeTrajectory(trajectory);
	int s;
	int d;

	adjustSpatialResolution(fieldOfViewXY, &trajectory->imageDimensions[0], &xySpatialResolution);
	trajectory->imageDimensions[1] = trajectory->imageDimensions[0];

	adjustSpatialResolution(fieldOfViewZ, &trajectory->imageDimensions[2], &zSpatialResolution);

	trajectory->spatialResolution[0] = xySpatialResolution;
	trajectory->spatialResolution[1] = trajectory->spatialResolution[0];
	trajectory->spatialResolution[2] = zSpatialResolution;
	trajectory->fieldOfView[0] = fieldOfViewXY;
	trajectory->fieldOfView[1] = trajectory->fieldOfView[0];
	trajectory->fieldOfView[2] = fieldOfViewZ;
	if(filterFieldOfView>0)
	{
		trajectory->maxReadoutGradientAmplitude = fminf(calculateMaxReadoutGradientAmplitude(filterFieldOfView, samplingInterval), maxGradientAmplitude);
	} 
	else 
		trajectory->maxReadoutGradientAmplitude = maxGradientAmplitude;
	trajectory->readoutPoints = ceil(readoutDuration/samplingInterval);
	trajectory->readoutPoints += trajectory->readoutPoints%2;
	trajectory->maxGradientAmplitude = maxGradientAmplitude;
	trajectory->maxSlewRate = maxSlewRate;
	trajectory->samplingInterval = samplingInterval;
	trajectory->dimensions = 3;
	trajectory->bases = bases;

	cones->rotatable = rotatable;
	cones->interconeCompensation = interConeCompensation;

	if(variableDensity)
	{
		trajectory->variableDensity = (struct VariableDensity*)malloc(sizeof(struct VariableDensity));
		copyVariableDensity(variableDensity, trajectory->variableDensity);
	}
	else
		trajectory->variableDensity = NULL;
	generateConesBasis(cones);
	//saveGradientWaveforms("grad.wav", cones->basisGradientWaveforms, 3, trajectory->bases, trajectory->waveformPoints, trajectory->readoutPoints, fmax(fieldOfViewXY,fieldOfViewZ), maxGradientAmplitude, 4, samplingInterval, "cones", LittleEndian);
	makeConesInterpolation(cones);
	trajectory->storage = storage;
	if(storage==StoreBasis)
	{
		trajectory->gradientWaveforms = (float*)malloc(bases*3*trajectory->waveformPoints*sizeof(float));
		trajectory->kSpaceCoordinates = (float*)malloc(bases*3*trajectory->readoutPoints*sizeof(float));
		memcpy(trajectory->gradientWaveforms, cones->basisGradientWaveforms, bases*3*trajectory->waveformPoints*sizeof(float));
		memcpy(trajectory->kSpaceCoordinates, cones->basisKspaceCoordinates, bases*3*trajectory->readoutPoints*sizeof(float));
	}
	else
	{
		trajectory->gradientWaveforms = (float*)malloc(trajectory->readouts*3*trajectory->waveformPoints*sizeof(float));
		trajectory->kSpaceCoordinates = (float*)malloc(trajectory->readouts*3*trajectory->readoutPoints*sizeof(float));
		for(s=0; s<cones->interpolation.readouts; s++)
		{
			generateReadoutWaveforms(s, cones, &trajectory->gradientWaveforms[s*3*trajectory->waveformPoints], &trajectory->gradientWaveforms[(s*3+1)*trajectory->waveformPoints], &trajectory->gradientWaveforms[(s*3+2)*trajectory->waveformPoints]);
			for(d=0; d<3; d++)
				gradientToKspace(&trajectory->gradientWaveforms[(s*3+d)*trajectory->waveformPoints], &trajectory->kSpaceCoordinates[(s*3+d)*trajectory->readoutPoints], trajectory->samplingInterval, trajectory->readoutPoints);
		}
	}

	return cones;
}

struct Cones *allocateCones(int bases)
{
	struct Cones *cones = (struct Cones*)malloc(sizeof(struct Cones));

	cones->coneAngles = NULL;
	cones->coneAngleDensityCompensation = NULL;
	cones->designConeAngles = NULL;
	cones->basisGradientWaveforms = NULL;
	cones->basisKspaceCoordinates = NULL;

	cones->basisReadoutPoints = (int*)malloc(bases*sizeof(int));
	cones->basisWaveformPoints = (int*)malloc(bases*sizeof(int));

	cones->trajectory.bases = bases;

	return cones;
}

void freeCones(struct Cones *cones)
{
	if(cones->coneAngles)
		free(cones->coneAngles);
	if(cones->coneAngleDensityCompensation)
		free(cones->coneAngleDensityCompensation);
	if(cones->designConeAngles)
		free(cones->designConeAngles);
	free(cones->basisReadoutPoints);
	free(cones->basisWaveformPoints);
	free(cones->basisGradientWaveforms);
	free(cones->basisKspaceCoordinates);
}

