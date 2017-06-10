/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University. 
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#include "angles.h"

#include "arrayops.h"
#include "mathops.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int MAX_ANGLES = 1000;

/*!
 * \brief Ellipse with diameters of X and Y
 * Peder Larson and Paul Gurney, 6/28/2006
 * \param angle	Angle at which to calculate extent
 * \param axisLengths	Ellipse axis lengths
 * \return diameter for input polar coordinates.
 */
float ellipseExtent(float angle, float *axisLengths)
{
	float X = axisLengths[0];
	float Y = axisLengths[1];

	float cosa = cos(angle);
	float sina  = sin(angle);

	return 1.0f/sqrt(cosa*cosa/(X*X) + sina*sina/(Y*Y));
}

float inverseEllipseExtent(float angle, float *axisLengths)
{
	float inverseAxisLengths[2];
	inverseAxisLengths[0] = 1/axisLengths[0];
	inverseAxisLengths[1] = 1/axisLengths[1];

	return 1/ellipseExtent(angle, inverseAxisLengths);
}

float getExtent(enum AngleShape shape, float angle, float *axisLengths)
{
	float extent=0;

	switch(shape)
	{
		case ConstantShape:
			extent = *axisLengths;
			break;
		case EllipticalShape:
			extent = ellipseExtent(angle, axisLengths);
			break;
		case InverseEllipticalShape:
			extent = inverseEllipseExtent(angle, axisLengths);
			break;
		default:
			fprintf(stderr, "Invalid shape function %d\n", shape);
			break;
	}

	return extent;
}

void calculateAngles(float theta0, float theta_width, enum AngleShape FOV, float *F, enum AngleShape KFCN, float *K, float **theta, float **kmax, float **dcf_theta, int *ntheta)
{
	int n;
	int N;
	float theta_cutoff;
	float dtheta_approx;
	float kmax_mid;
	float dtheta;

	float thetaTemp[MAX_ANGLES];
	float kmaxTemp[MAX_ANGLES];
	float dcfTemp[MAX_ANGLES];

	int allocbuf = 2;

n = 0;
thetaTemp[n] = theta0;
theta_cutoff = theta0 + theta_width;

while (thetaTemp[n] < theta_cutoff)
{
  /*dtheta_approx = 1 / ( feval(KFCN, theta(n), K{:})/2 * feval(FOV, theta(n) + pi/2, F{:}));*/
	dtheta_approx = 1.0f / ( getExtent(KFCN, thetaTemp[n], K)/2.0f * getExtent(FOV, thetaTemp[n] + M_PI_2, F));
  /*kmax_mid = feval(KFCN, theta(n) + dtheta_approx/2, K{:})/2;*/
	kmax_mid = getExtent(KFCN, thetaTemp[n] + dtheta_approx/2.0f, K)/2.0f;

 /* dtheta = 1 / ( kmax_mid * feval(FOV, theta(n) + dtheta_approx/2 + pi/2, F{:}));*/
	dtheta = 1 / ( kmax_mid * getExtent(FOV, thetaTemp[n] + dtheta_approx/2 + M_PI_2, F));
  thetaTemp[n+1] = thetaTemp[n] + dtheta;
  n++;
}

N = n+1;

/*% adjust theta for symmetry
% choose adjustment based on which spoke is closest to pi
% NOTE: could also choose adjustment based on whether an even or odd number
% of full-spokes is desired*/
if ((thetaTemp[N-1] - theta_cutoff) > (theta_cutoff - thetaTemp[N-2]))
{
	N -= 2;
  /*theta = (theta(1:end-2) - theta0)*theta_width/(theta(end-1) - theta0) + theta0;*/

}
else
{
	N--;
  /*theta = (theta(1:end-1) - theta0)*theta_width/(theta(end) - theta0) + theta0;*/
}

addfloat(thetaTemp, -theta0, N);
scalefloats(thetaTemp, N, theta_width/(thetaTemp[N]-theta0));
addfloat(thetaTemp, theta0, N);

/* Calculate kmax function
kmax = feval(KFCN, theta, K{:})/2;*/

/* DCF = projection length * projection spacing
dcf_theta = kmax./feval(FOV, theta + pi/2, F{:});*/
for(n=0; n<N; n++)
{
	kmaxTemp[n] = getExtent(KFCN, thetaTemp[n], K)/2;

	dcfTemp[n] = kmaxTemp[n]/getExtent(FOV, thetaTemp[n] + M_PI_2, F);
}

	if(theta)
	{
		*theta = (float*)malloc((N+1+allocbuf)*sizeof(float));
		memcpy(*theta, thetaTemp, N*sizeof(float));
	}
	if(kmax)
	{
		*kmax = (float*)malloc((N+1+allocbuf)*sizeof(float));
		memcpy(*kmax, kmaxTemp, N*sizeof(float));
	}
	if(dcf_theta)
	{
		*dcf_theta = (float*)malloc((N+1+allocbuf)*sizeof(float));
		memcpy(*dcf_theta, dcfTemp, N*sizeof(float));
	}

	*ntheta = N;

	/*for(n=0; n<2; n++)
		K[n] *= .5;*/

	return;
}

void calculateAngles3D(int proj_type, enum AngleShape thetaShape, enum AngleShape phiShape, int *N, float **thetaOut, float **phiOut, float **kmaxOut, float **dcfOut, int *Nproj)
{
	float Ft[2];	/* yz FOV */
	float Fp[2];	/* transverse FOV */
	float Kt[2];
	float thwid;
	/*float (*FOV)(float,float*);
	float (*FOVtheta)(float,float*);
	float (*FOVphi)(float,float*);*/
//	float (*KFCNtheta)(float,float*);
	float *theta;
	float *phi;
	float *kmax;
	float *dcf;
	/*float thcones[MAXANG];
	float kmaxcones[MAXANG];
	float dcones[MAXANG];*/
	float *thcones;
	float *kmaxcones;
	/*float *dcones;*/
	int k;
	float rk;
	float dphi_approx;
	float dphi;
	float *t;
	float *n;
	float kmaxest;
	int ncones;
	int Nphiest;
	int Np;

	enum AngleShape kShape = ConstantShape;


/* FOVphi */
	Ft[0] = N[2];
	Ft[1] = N[0];

//	KFCNtheta = &shapeConst;
	Kt[0] = 1.0f;

	Fp[0] = N[0];
	Fp[1] = N[1];

/* Wong & Roos style 3D PR design */

if(proj_type) /* half projections */
  thwid = M_PI;
else /* full projections */
  thwid = M_PI_2;

/*[thcones, kmaxcones, dcones] = calc_angles(0, thwid, FOVtheta, Ft{:}, KFCNtheta, Kt{:});*/
calculateAngles(0.0f, thwid, thetaShape, Ft, kShape, Kt, &thcones, &kmaxcones, NULL, &ncones);

thcones[ncones] = thwid;
/*kmaxcones[Ncones] = feval(KFCNtheta, thwid, Kt{:})/2;*/
kmaxcones[ncones] = getExtent(ConstantShape, thwid, Kt)/2;
ncones++;

kmaxest = getExtent(ConstantShape, M_PI_2, Kt);
/*phiest = calc_angles(0, 2*pi, FOVphi, Fp{:}, @const, kmaxest);*/
calculateAngles(0, 2*M_PI, phiShape, Fp, kShape, &kmaxest, NULL, NULL, NULL, &Nphiest);
/*Nphiest = length(phiest);*/

/*ncones = length(thcones);*/

/*for k = 2:ncones*/
t = (float*)malloc((ncones+1)*sizeof(float));
t[0] = 1;
for(k=1; k<ncones; k++)
  t[k] = t[k-1] + Nphiest * sin((thcones[k] + thcones[k-1])/2) * (kmaxcones[k] + kmaxcones[k-1]) / kmaxest;

if (proj_type == 0)
{
  /* add extra quarter turn of spiral */
  thcones[ncones] = M_PI_2 + 1/(kmaxest/2 * getExtent(thetaShape, M_PI, Ft)) / 4.0f;

  t[ncones] = t[ncones-1] + Nphiest/4.0f;

  kmaxcones[ncones] = getExtent(ConstantShape, thcones[ncones-1], Kt)/2.0f;
}

ncones++;

/*n = 1:floor(t(end));*/
Np = t[ncones-1];
n = (float*)malloc(Np*sizeof(float));
*thetaOut = (float*)malloc(Np*sizeof(float));
*kmaxOut = (float*)malloc(Np*sizeof(float));
*phiOut = (float*)malloc(Np*sizeof(float));
*dcfOut = (float*)malloc(Np*sizeof(float));

for(k=0; k<Np; k++)
	n[k] = k+1;

theta = *thetaOut;
/*theta = interp1(t, thcones, n, 'linear');*/
interpolatefloats(t, thcones, ncones, n, theta, Np);

kmax = *kmaxOut;
/*kmax = interp1(t, kmaxcones, n, 'linear');*/
interpolatefloats(t, kmaxcones, ncones, n, kmax, Np);

phi = *phiOut;
phi[0] = 0;
/*for k = 2:length(n)*/
for(k=1; k<Np; k++)
{
  rk = kmax[k] * sin(theta[k]);

   /*dphi_approx = 1 / (rk * feval(FOVphi, phi(k-1) + pi/2, Fp{:}));*/
		dphi_approx = 1.0f/(rk*getExtent(phiShape, phi[k-1]+M_PI_2, Fp));

  /*dphi = 1 / (rk * feval(FOVphi, phi(k-1) + dphi_approx/2  + pi/2, Fp{:}));*/
		dphi = 1.0f / (rk * getExtent(phiShape, phi[k-1] + dphi_approx/2.0f  + M_PI_2, Fp));

  phi[k] = phi[k-1] + dphi;
}

dcf = *dcfOut;
for(k=0; k<Np; k++)
{
	/*dcf = (kmax./feval(FOVtheta, theta + pi/2, Ft{:})) ./ feval(FOVphi, phi + pi/2, Fp{:});*/
	dcf[k] = (kmax[k]/getExtent(thetaShape, theta[k] + M_PI_2, Ft)) / getExtent(phiShape, phi[k] + M_PI_2, Fp);
}

if (proj_type == 0) /* adjust dcf for full projections near k_xy */
{
	for(k=0; k<Np; k++)
	{
  /*I1 = find((phi > phi(end) - 2*pi) & (phi <= phi(end) - pi));
  dcf(I1) = dcf(I1) .* (1 - .5 * (phi(I1)- (phi(end) - 2*pi)) / pi );*/
		if((phi[k]>phi[Np-1]-2.0*M_PI) && (phi[k] <= phi[Np-1]-M_PI))
			dcf[k] = dcf[k] * (1 - 0.5f * (phi[k] - (phi[Np-1] - 2.0*M_PI)) * M_1_PI );
  /*I2 = find(phi > phi(end) - pi);
  dcf(I2) = dcf(I2) .* (1 - .5 * (phi(I2)- (phi(end) - pi)) / pi);*/
		if(phi[k] > (phi[Np-1]-M_PI))
			dcf[k] = dcf[k] * (1 - .5f * (phi[k]- (phi[Np-1] - M_PI)) * M_1_PI);

		phi[k] = modulus(phi[k], 2.0*M_PI);
	}
}

/*phi = mod(phi, 2*pi);*/
	*Nproj = Np;

	free(n);

	return;
}
