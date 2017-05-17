/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University. 
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

**************************************************************************/
#include "angles.h"

#include "arrayops.h"

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

/*% Make theta0, theta_width optional
if isa(theta0, 'function_handle')
  varargin(3:end+2) = varargin;
  varargin{1} = theta_width;
  if nargin > 2
	varargin{2} = FOV;
  end
  FOV = theta0;
  theta0 = 0;
  theta_width = pi;
elseif isa(theta_width, 'function_handle')
  varargin(2:end+1) = varargin;
  varargin{1} = FOV;
  FOV = theta_width;
  theta_width = pi;
end*/

/*for k = 1:length(varargin)
  if isa(varargin{k}, 'function_handle')
	KFCN = varargin{k};
	K = varargin(k+1:end);
	break
  else
	F(k) = varargin(k);
  end
end

if (~exist('KFCN'))
  KFCN = @const;
  K = {1};
end*/

	/*for(n=0; n<2; n++)
		K[n] *= 2;*/

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

