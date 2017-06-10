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

void grd_triangle(float GMAX, float SMAX, float Ktgt, float Ts, float** kx, float**gx, float**slew, float**time, int *N, int *N_ramp);

void grd_trapezoid(float GMAX, float SMAX, float Ktgt, float Ts, float** kx, float** gx, float** slew, float** time, int *N, int *Nrout);

/**
% Project: Concentric Rings
% filename: grd_rampup.m
%           rampup gradient to target strength
% Modified:
%           2009/02/28 copied from rampup_func.m
% 2009/02/28
% Hochong Wu

function [kx,gx,slew,time] = grd_rampup( GYROMAGNETIC_RATIO,GMAX,SMAX, Gtgt,Ts, DEBUG_MRGRADIENT)
% % Hardware constraints
% TMIN = 4 *10^-6;            % 4 us --> 4e-6 s
% GMAX = 40 *(10/100);        % 40 mT/m --> 4 G/cm
% SMAX = 150 *(10/100/10^-3); % 150 mT/m/ms --> 15,000 G/cm/s
% % Constants
% GYROMAGNETIC_RATIO = 4258; % Hz/G
*/
void grd_rampup(float GMAX, float SMAX, float Gtgt, float Ts, float** kx, float** gx, float** slew, float** time, int *N)
{
	float T_ramp;
	int N_ramp;
	float scaling;

	int n;

/* Check */
if( Gtgt > 1.001*GMAX )
	printf("grd_rampup: Gtgt > 1.001*GMAX\n");

/* Ramp timing */
T_ramp = Gtgt/SMAX; /* sec */
/* convert to number of samples
% use 1 extra sample(s), can scale down at the end*/
N_ramp = ceil( T_ramp/Ts ) + 1;

*gx = (float*)malloc(N_ramp*sizeof(float));

/* Ramp up
gx = [0:N_ramp-1]*Ts*SMAX;*/
for(n=0; n<N_ramp; n++)
	(*gx)[n] = n*Ts*SMAX;

/* Scale the ramp to have exactly GMAX at the end */
scaling = (Gtgt/(*gx)[N_ramp-1]);
if(fabs(scaling)>1)
	printf("grd_rampup: abs(scaling)=%.2f> 1 !\n", scaling);

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
% filename: grd_triangle.m
%           generate triangle gradient lobes
% Modified:
%           2009/02/28 copied from triangle_func.m
% 2009/02/28
% Hochong Wu

function [kx, gx, slew, time] = grd_triangle( GYROMAGNETIC_RATIO,GMAX,SMAX, Ktgt,Ts, DEBUG_MRGRADIENT )*/
void grd_triangle(float GMAX, float SMAX, float Ktgt, float Ts, float** kx, float**gx, float**slew, float**time, int *N, int *N_ramp)
{
/*% % Hardware constraints
% TMIN = 4 *10^-6;            % 4 us --> 4e-6 s
% GMAX = 40 *(10/100);        % 40 mT/m --> 4 G/cm
% SMAX = 150 *(10/100/10^-3); % 150 mT/m/ms --> 15,000 G/cm/s
% % Constants
% GYROMAGNETIC_RATIO = 4258; % Hz/G*/

	float T_ramp;
	/*int N_ramp;*/
//	int N_const = 0;
	int Ngx;
	float scaling;

	int n;

	/* Check */
	if( fabs(Ktgt/GYROMAGNETIC_RATIO) > (GMAX*GMAX)/SMAX )
	{
		if( DEBUG_MRGRADIENT)
		   printf("grd_triangle: using a trapezoid lobe\n");
		grd_trapezoid(GMAX, SMAX, Ktgt,Ts, kx, gx, slew, time, N, N_ramp);
		return;
	}

	/* Ramp timing */
	T_ramp = sqrt( fabs(Ktgt/GYROMAGNETIC_RATIO)/SMAX ); /* sec  */
	/*% convert to number of samples
	% use 3 extra sample(s), can scale down at the end*/
	*N_ramp = ceil( T_ramp/Ts ) + 3;

	/* The max amplitude of the triangle, not explicitly used *
	%Gmax = abs(Ktgt)/GYROMAGNETIC_RATIO/T_ramp;

	% 1. Ramp up
	rampup = [0:N_ramp-1]*Ts*SMAX;*/
	Ngx = 2**N_ramp;
	*gx = (float*)malloc(Ngx*sizeof(float));

	for(n=0; n<*N_ramp; n++)
	{
		(*gx)[n] = n*Ts*SMAX;
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
		printf("grd_triangle: abs(scaling)=%.2f> 1 !\n", scaling);

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
% filename: grd_trapezoid.m
%           generate trapezoid gradient lobes
% Modified:
%           2009/02/28 copied from trapezoid_func.m
% 2009/02/28
% Hochong Wu

function [kx, gx, slew, time] = grd_trapezoid( GYROMAGNETIC_RATIO,GMAX,SMAX, Ktgt,Ts, DEBUG_MRGRADIENT )
% % Hardware constraints
% TMIN = 4 *10^-6;            % 4 us --> 4e-6 s
% GMAX = 40 *(10/100);        % 40 mT/m --> 4 G/cm
% SMAX = 150 *(10/100/10^-3); % 150 mT/m/ms --> 15,000 G/cm/s
% % Constants
% GYROMAGNETIC_RATIO = 4258; % Hz/G
*/

void grd_trapezoid(float GMAX, float SMAX, float Ktgt, float Ts, float** kx, float** gx, float** slew, float** time, int *N, int *Nrout)
{

	float T_ramp;
	float T_plateau;
	int Nr;
	int Npl;
	int Ngx;
	float scaling;

	int n;

/* Check */
if( fabs(Ktgt) < GYROMAGNETIC_RATIO*(GMAX*GMAX)/SMAX )
{
	grd_triangle(GMAX,SMAX, Ktgt,Ts, kx, gx, slew, time, N, &n);
	return;
}

/* Ramp timing */
T_ramp = GMAX/SMAX; /* sec */
T_plateau = fabs(Ktgt)/GMAX/GYROMAGNETIC_RATIO - T_ramp; /* sec*/

/* 1. Ramp up
rampup = [0:Ts:T_ramp-Ts]*SMAX;*/
Nr = T_ramp/Ts;

/* 2. Plateau, use 3 extra sample(s), can scale down at the end
plateau = ones(1, ceil(T_plateau/Ts)+3) * rampup(end);*/
Npl = ceil(T_plateau/Ts)+3;

/* 3. Ramp down
rampdown = plateau(end) - rampup;

gx = [rampup plateau rampdown];
kx = GYROMAGNETIC_RATIO*Ts*cumsum(gx);*/

Ngx = 2*Nr+Npl;

*gx = (float*)malloc(Ngx*sizeof(float));

for(n=0; n<Nr; n++)
	(*gx)[n] = n*Ts*SMAX;
for(n=Nr; n<(Nr+Npl); n++)
	(*gx)[n] = (*gx)[Nr-1];
for(n=0; n<Nr; n++)
	(*gx)[Nr+Npl+n] = (*gx)[Nr-1-n];

/* Scale the trapezoid to have exactly Ktgt area */
scaling = Ktgt/(sumfloats((*gx), Ngx)*GYROMAGNETIC_RATIO*Ts);
if( fabs(scaling)>1 && DEBUG_MRGRADIENT )
	printf("grd_trapezoid: abs(scaling)=%.2f> 1 !\n", scaling);

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


void grd_readout(int doDephase, int doRephase, int nread, float res, float gmax, float smax, float Ts, float **g, int *nramp, int *ndep, int *npts)
{
	float G;
	float kmax = 5/res;

	float *gramp, *gtrap;
	float kramp;

	int ntrap;

	int i;
	int offset;

	G = 2*kmax/(GYROMAGNETIC_RATIO*nread*Ts);

	grd_rampup(gmax, smax, G, Ts, NULL, &gramp, NULL, NULL, nramp);
	kramp = GYROMAGNETIC_RATIO*Ts*sumfloats(gramp, *nramp);
	grd_trapezoid(gmax, smax, -kmax-kramp, Ts, NULL, &gtrap, NULL, NULL, &ntrap, &i);

	/*npts = (*ndep)*2 + nread;*/
	*npts = (doDephase>0)*(ntrap) + *nramp + nread + *nramp + (doRephase>0)*(ntrap);

	*g = (float*)malloc(*npts*sizeof(float));

	if(doDephase)
	{
//		copyfloats(gtrap, *g, ntrap);
		memcpy(*g, gtrap, ntrap*sizeof(float));
		offset = ntrap;
		*ndep = ntrap + *nramp;
	}
	else
	{
		*ndep = *nramp;
		offset = 0;
	}

	/* copy ramp up */
//	copyfloats(gramp, &(*g)[offset], *nramp);
	memcpy(&(*g)[offset], gramp, *nramp*sizeof(float));
	offset += *nramp;

	/* copy readout portion */
	for(i=0; i<nread; i++)
	   (*g)[i+offset] = G;
	offset += nread;

	/* copy ramp down */
	for(i=0; i<*nramp; i++)
		(*g)[i+offset] = gramp[*nramp-i-1];
	offset += *nramp;

	if(doRephase)
		/* copy rephaser */
		for(i=0; i<ntrap; i++)
		   /*(*g)[i+offset] = (*g)[*ndep-1-i];*/
			(*g)[i+offset] = gtrap[i];

	free(gramp);
	free(gtrap);

	return;
}
