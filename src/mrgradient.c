#include "mrgradient.h"

#ifndef GYROMAGNETIC_RATIO
#define GYROMAGNETIC_RATIO 4257
#endif

#include "arrayops.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define DEBUG_MRGRADIENT 1

void triangleGradient(float gradientLimit, float slewRateLimit, float Ktgt, float Ts, float** kx, float**gx, float**slew, float**time, int *N, int *N_ramp);

void trapezoidGradient(float gradientLimit, float slewRateLimit, float Ktgt, float Ts, float** kx, float** gx, float** slew, float** time, int *N, int *Nrout);

/**
% Project: Concentric Rings
% filename: rampUpGradient.m
%           rampup gradient to target strength
% Modified:
%           2009/02/28 copied from rampup_func.m
% 2009/02/28
% Hochong Wu

function [kx,gx,slew,time] = rampUpGradient( GYROMAGNETIC_RATIO,gradientLimit,slewRateLimit, gradientTarget,Ts, DEBUG_MRGRADIENT)
% % Hardware constraints
% TMIN = 4 *10^-6;            % 4 us --> 4e-6 s
% gradientLimit = 40 *(10/100);        % 40 mT/m --> 4 G/cm
% slewRateLimit = 150 *(10/100/10^-3); % 150 mT/m/ms --> 15,000 G/cm/s
% % Constants
% GYROMAGNETIC_RATIO = 4258; % Hz/G
*/
void rampGradient(float slewRateLimit, float gradientInitial, float gradientTarget, float Ts, float** kx, float** gx, float** slew, float** time, int *N)
{
	float duration;
	int N_ramp;
	float scaling;

	int n;

	/* Ramp timing */
	duration = fabs(gradientTarget-gradientInitial)/slewRateLimit; /* sec */
	/* convert to number of samples
	% use 1 extra sample(s), can scale down at the end*/
	N_ramp = ceil( duration/Ts ) + 1;

	*gx = (float*)malloc(N_ramp*sizeof(float));

	/* Ramp up
	gx = [0:N_ramp-1]*Ts*slewRateLimit;*/
	for(n=0; n<N_ramp; n++)
		(*gx)[n] = n*Ts*slewRateLimit;

	/* Scale the ramp to have exactly gradientLimit at the end */
	scaling = (gradientTarget/(*gx)[N_ramp-1]);
	if(fabs(scaling)>1)
		printf("rampUpGradient: abs(scaling)=%.2f> 1 !\n", scaling);

	/* Update output
	gx   = gx * scaling;*/
	scalefloats(*gx, N_ramp, scaling);

	/*kx   = GYROMAGNETIC_RATIO*Ts*cumsum(gx);*/
	if(kx)
	{
		*kx = (float*)malloc(N_ramp*sizeof(float));
		cumsum(*gx, *kx, N_ramp);
		scalefloats(*kx, N_ramp, GYROMAGNETIC_RATIO*Ts);
	}

	/*slew = diff([0 gx])/Ts;
	time = Ts*[0:numel(gx)-1];*/
	if(slew)
	{
		*slew = (float*)malloc(N_ramp*sizeof(float));
		diffArray(*gx, *slew, N_ramp);
		scalefloats(*slew, N_ramp, 1.0f/Ts);
	}

	if(time)
	{
		*time = (float*)malloc(N_ramp*sizeof(float));
		for(n=0; n<N_ramp; n++)
			(*time)[n] = n*Ts;
	}

	*N = N_ramp;

	return;
}

/**
% Project: Concentric Rings
% filename: triangleGradient.m
%           generate triangle gradient lobes
% Modified:
%           2009/02/28 copied from triangle_func.m
% 2009/02/28
% Hochong Wu

function [kx, gx, slew, time] = triangleGradient( GYROMAGNETIC_RATIO,gradientLimit,slewRateLimit, Ktgt,Ts, DEBUG_MRGRADIENT )*/
void triangleGradient(float gradientLimit, float slewRateLimit, float Ktgt, float Ts, float** kx, float**gx, float**slew, float**time, int *N, int *N_ramp)
{
/*% % Hardware constraints
% TMIN = 4 *10^-6;            % 4 us --> 4e-6 s
% gradientLimit = 40 *(10/100);        % 40 mT/m --> 4 G/cm
% slewRateLimit = 150 *(10/100/10^-3); % 150 mT/m/ms --> 15,000 G/cm/s
% % Constants
% GYROMAGNETIC_RATIO = 4258; % Hz/G*/

	float duration;
	/*int N_ramp;*/
	int Ngx;
	float scaling;

	int n;

	/* Check */
	if( fabs(Ktgt/GYROMAGNETIC_RATIO) > (gradientLimit*gradientLimit)/slewRateLimit )
	{
		if( DEBUG_MRGRADIENT)
		   printf("triangleGradient: using a trapezoid lobe\n");
		trapezoidGradient(gradientLimit, slewRateLimit, Ktgt,Ts, kx, gx, slew, time, N, N_ramp);
		return;
	}

	/* Ramp timing */
	duration = sqrt( fabs(Ktgt/GYROMAGNETIC_RATIO)/slewRateLimit ); /* sec  */
	/*% convert to number of samples
	% use 3 extra sample(s), can scale down at the end*/
	*N_ramp = ceil( duration/Ts ) + 3;

	/* The max amplitude of the triangle, not explicitly used *
	%gradientLimit = abs(Ktgt)/GYROMAGNETIC_RATIO/duration;

	% 1. Ramp up
	rampup = [0:N_ramp-1]*Ts*slewRateLimit;*/
	Ngx = 2**N_ramp;
	*gx = (float*)malloc(Ngx*sizeof(float));

	for(n=0; n<*N_ramp; n++)
	{
		(*gx)[n] = n*Ts*slewRateLimit;
		(*gx)[Ngx-1-n] = (*gx)[n];
	}

	/*% 2. Ramp down
	rampdown = fliplr( rampup(1:end-1) );*/

	/*gx = [rampup rampdown];
	kx = GYROMAGNETIC_RATIO*Ts*cumsum(gx);*/

	/* Scale the triangle to have exactly Ktgt area*/
	/*scaling = (Ktgt/(*kx)[Ngx-1]);*/
	scaling = Ktgt/(sumfloats((*gx), Ngx)*GYROMAGNETIC_RATIO*Ts);
	if( fabs(scaling)>1 && DEBUG_MRGRADIENT )
		printf("triangleGradient: abs(scaling)=%.2f> 1 !\n", scaling);

	/* Update output */
	/*gx   = gx * scaling;*/
	scalefloats(*gx, Ngx, scaling);

	/*kx   = kx * scaling;*/
	if(kx)
	{
		*kx = (float*)malloc(Ngx*sizeof(float));
		cumsum(*gx, *kx, Ngx);
		scalefloats(*kx, Ngx, GYROMAGNETIC_RATIO*Ts);
	}

	/*slew = diff([0 gx])/Ts;*/
	if(slew)
	{
		*slew = (float*)malloc(Ngx*sizeof(float));
		diffArray(*gx, *slew, Ngx);
		scalefloats(*slew, Ngx, 1/Ts);
	}

	/*time = Ts*[0:numel(gx)-1];*/
	if(time)
	{
		*time = (float*)malloc(Ngx*sizeof(float));
	for(n=0; n<Ngx; n++)
		(*time)[n] = n*Ts;
	}

	*N = Ngx;

	return;
}

/**
% Project: Concentric Rings
% filename: trapezoidGradient.m
%           generate trapezoid gradient lobes
% Modified:
%           2009/02/28 copied from trapezoid_func.m
% 2009/02/28
% Hochong Wu

function [kx, gx, slew, time] = trapezoidGradient( GYROMAGNETIC_RATIO,gradientLimit,slewRateLimit, Ktgt,Ts, DEBUG_MRGRADIENT )
% % Hardware constraints
% TMIN = 4 *10^-6;            % 4 us --> 4e-6 s
% gradientLimit = 40 *(10/100);        % 40 mT/m --> 4 G/cm
% slewRateLimit = 150 *(10/100/10^-3); % 150 mT/m/ms --> 15,000 G/cm/s
% % Constants
% GYROMAGNETIC_RATIO = 4258; % Hz/G
*/

void trapezoidGradient(float gradientLimit, float slewRateLimit, float Ktgt, float Ts, float** kx, float** gx, float** slew, float** time, int *N, int *Nrout)
{

	float duration;
	float durationPlateau;
	int Nr;
	int Npl;
	int Ngx;
	float scaling;

	int n;
	if( fabs(Ktgt) < GYROMAGNETIC_RATIO*(gradientLimit*gradientLimit)/slewRateLimit )
	{
		triangleGradient(gradientLimit,slewRateLimit, Ktgt,Ts, kx, gx, slew, time, N, &n);
		return;
	}

	duration = gradientLimit/slewRateLimit;
	durationPlateau = fabs(Ktgt)/gradientLimit/GYROMAGNETIC_RATIO - duration;

	/* 1. Ramp up
	rampup = [0:Ts:duration-Ts]*slewRateLimit;*/
	Nr = duration/Ts;

	/* 2. Plateau, use 3 extra sample(s), can scale down at the end
	plateau = ones(1, ceil(durationPlateau/Ts)+3) * rampup(end);*/
	Npl = ceil(durationPlateau/Ts)+3;

	/* 3. Ramp down
	rampdown = plateau(end) - rampup;

	gx = [rampup plateau rampdown];
	kx = GYROMAGNETIC_RATIO*Ts*cumsum(gx);*/

	Ngx = 2*Nr+Npl;

	*gx = (float*)malloc(Ngx*sizeof(float));

	for(n=0; n<Nr; n++)
		(*gx)[n] = n*Ts*slewRateLimit;
	for(n=Nr; n<(Nr+Npl); n++)
		(*gx)[n] = (*gx)[Nr-1];
	for(n=0; n<Nr; n++)
		(*gx)[Nr+Npl+n] = (*gx)[Nr-1-n];

	/* Scale the trapezoid to have exactly Ktgt area */
	scaling = Ktgt/(sumfloats((*gx), Ngx)*GYROMAGNETIC_RATIO*Ts);
	if( fabs(scaling)>1 && DEBUG_MRGRADIENT )
		printf("trapezoidGradient: abs(scaling)=%.2f> 1 !\n", scaling);

	/* Update output */
	/*gx   = gx * scaling;
	kx   = kx * scaling;
	slew = diff([0 gx])/Ts;
	time = Ts*[0:numel(gx)-1];*/

	scalefloats(*gx, Ngx, scaling);

	if(kx)
	{
		*kx = (float*)malloc(Ngx*sizeof(float));
		cumsum(*gx, *kx, Ngx);
		scalefloats(*kx, Ngx, GYROMAGNETIC_RATIO*Ts);
	}

	if(slew)
	{
		*slew = (float*)malloc(Ngx*sizeof(float));
		diffArray(*gx, *slew, Ngx);
		scalefloats(*slew, Ngx, 1/Ts);
	}

	if(time)
	{
		*time = (float*)malloc(Ngx*sizeof(float));
		for(n=0; n<Ngx; n++)
			(*time)[n] = n*Ts;
	}

	*N = Ngx;
	*Nrout = Nr;

	return;
}


void readoutGradient(int doDephase, int doRephase, int nread, float spatialResolution, float gradientLimit, float slewRateLimit, float Ts, float **g, int *rampPoints, int *ndep, int *npts)
{
	float kmax = 5/spatialResolution;

	float *gramp, *gtrap;
	float kramp;

	int ntrap;

	int i;
	int offset;

	float readoutGradientAmplitude = 2*kmax/(GYROMAGNETIC_RATIO*nread*Ts);

	rampGradient(slewRateLimit, 0, readoutGradientAmplitude, Ts, NULL, &gramp, NULL, NULL, rampPoints);
	kramp = GYROMAGNETIC_RATIO*Ts*sumfloats(gramp, *rampPoints);
	trapezoidGradient(gradientLimit, slewRateLimit, -kmax-kramp, Ts, NULL, &gtrap, NULL, NULL, &ntrap, &i);

	/*npts = (*ndep)*2 + nread;*/
	*npts = (doDephase>0)*(ntrap) + *rampPoints + nread + *rampPoints + (doRephase>0)*(ntrap);

	*g = (float*)malloc(*npts*sizeof(float));

	if(doDephase)
	{
		memcpy(*g, gtrap, ntrap*sizeof(float));
		offset = ntrap;
		*ndep = ntrap + *rampPoints;
	}
	else
	{
		*ndep = *rampPoints;
		offset = 0;
	}

	/* copy ramp up */
	memcpy(&(*g)[offset], gramp, *rampPoints*sizeof(float));
	offset += *rampPoints;

	/* copy readout portion */
	for(i=0; i<nread; i++)
	   (*g)[i+offset] = readoutGradientAmplitude;
	offset += nread;

	/* copy ramp down */
	for(i=0; i<*rampPoints; i++)
		(*g)[i+offset] = gramp[*rampPoints-i-1];
	offset += *rampPoints;

	if(doRephase)
		/* copy rephaser */
		for(i=0; i<ntrap; i++)
		   /*(*g)[i+offset] = (*g)[*ndep-1-i];*/
			(*g)[i+offset] = gtrap[i];

	free(gramp);
	free(gtrap);

	return;
}

void spoilerGradient(float gradientLimit, float slewRateLimit, float gradientInitial, float deltaKspace, float gradientFinal, float samplingInterval, float** gradientWaveform, int* points)
{
	int useTriangle = 1;

	float gradientTriangleAmplitude = slewRateLimit*deltaKspace/GYROMAGNETIC_RATIO+.5*(gradientInitial*gradientInitial-gradientFinal*gradientFinal);
	if(gradientTriangleAmplitude<0)
		useTriangle = 0;
	else
	{
		gradientTriangleAmplitude = sqrt(gradientTriangleAmplitude);
		if(gradientTriangleAmplitude>gradientLimit)
			useTriangle = 0;
	}

	if(useTriangle)
		gradientLimit = fmin(gradientTriangleAmplitude, gradientLimit);

	float areaRampUp = (gradientLimit*gradientLimit - gradientInitial*gradientInitial)/(2*slewRateLimit)*GYROMAGNETIC_RATIO;

	float areaRampDown;
	float areaRewind = 0;
	if(gradientFinal>=0)
	{
		areaRampDown = (gradientLimit*gradientLimit - gradientInitial*gradientInitial)/(2*slewRateLimit);
	}
	else
	{
		areaRampDown = gradientLimit*gradientLimit/(2*slewRateLimit)*GYROMAGNETIC_RATIO;
		areaRewind = gradientFinal*gradientFinal/(2*slewRateLimit)*GYROMAGNETIC_RATIO;
	}

	float areaPlateau = deltaKspace - areaRampUp - areaRampDown + areaRewind;
	float plateauTime = areaPlateau/(gradientLimit*GYROMAGNETIC_RATIO);

	float* gradientWaveformRampUp;
	int pointsRampUp;

	rampGradient(slewRateLimit, gradientInitial, gradientLimit, samplingInterval, NULL, &gradientWaveformRampUp, NULL, NULL, &pointsRampUp);

	int pointsPlateau;
	float* gradientWaveformPlateau;
	if(useTriangle)
	{
		pointsPlateau = 0;
	}
	else
	{
		pointsPlateau = plateauTime/samplingInterval;
		gradientWaveformPlateau = (float*)malloc(pointsPlateau*sizeof(float));
		int n;
		for(n=0; n<pointsPlateau; n++)
			gradientWaveformPlateau[n] = gradientLimit;
	}

	float* gradientWaveformRampDown;
	int pointsRampDown;

	float rampDownGradientAmplitude;
	if(gradientInitial*gradientFinal>0)
		rampDownGradientAmplitude = gradientFinal;
	else
		rampDownGradientAmplitude = 0;

	rampGradient(slewRateLimit, gradientLimit, rampDownGradientAmplitude, samplingInterval, NULL, &gradientWaveformRampDown, NULL, NULL, &pointsRampDown);

	float* gradientWaveformRewind;
	int pointsRewind = 0;
	if(gradientFinal)
		rampGradient(slewRateLimit, 0, gradientFinal, samplingInterval, NULL, &gradientWaveformRewind, NULL, NULL, &pointsRewind);

	*points = pointsRampUp + pointsPlateau + pointsRampDown + pointsRewind;
	*gradientWaveform = (float*)malloc(*points*sizeof(float));
	memcpy(*gradientWaveform, gradientWaveformRampUp, pointsRampUp*sizeof(float));
	if(!useTriangle)
	{
		memcpy(&(*gradientWaveform)[pointsRampDown], gradientWaveformPlateau, pointsPlateau*sizeof(float));
	}
	memcpy(&(*gradientWaveform)[pointsRampDown+pointsPlateau], gradientWaveformRampUp, pointsRampDown*sizeof(float));
	if(pointsRewind)
		memcpy(&(*gradientWaveform)[pointsRampDown+pointsPlateau+pointsRampDown], gradientWaveformRewind, pointsRewind*sizeof(float));
}
