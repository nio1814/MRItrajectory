#include "mathops.h"

#include <math.h>


int roundToInteger(float num)
{
	int out;
	int itg;
	float dec;

	itg = (int)num;
	dec = num-itg;
	if(dec<=0.5)
		out = (int)floor(num);
	else
		out = (int)ceil(num);

	return out;
}

void calculatePhase(float *realdata, float *imagdata, float *phasedata, long npts, float phasestart, float phasescale)
{
	long count;
	float radianphase;
	phasescale = 1.0;

	/*	printf("atan2 of (-2,1) = %g \n",atan2(-2.0,1.0)); */

	for (count = 0; count < npts; count++)
	{
		radianphase = atan2(*imagdata++, *realdata++);
		*phasedata++ = (radianphase+phasestart+3.14159265) * phasescale;		
	}

	return;
}

void calculateMagnitude(float *realdata, float *imagdata, float *magdata, long npts)
{
	long count;
	float magval;

	for (count=0; count < npts; count++)
	{
		magval = (*realdata) * (*realdata) + (*imagdata) * (*imagdata);
		*magdata++ = sqrt(magval);	

		realdata++;
		imagdata++;
	}
}

void unwrapPhase(float *phasedata, float* uPhaseData, int numPoints, int frsize, float tol)
{
	int i;
	float nextPhase, currentPhase;
	float phaseDiff;
	static float offset = 0.0f;

	offset = -phasedata[0];
	uPhaseData[0] = 0;

	for(i=0; i<(numPoints-1); i++)
	{
		nextPhase = phasedata[i+1];
		currentPhase = phasedata[i];
		phaseDiff = nextPhase-currentPhase;

		if((i+1)%frsize == 0)
			offset = -phasedata[i+1];
		else
		{
			if(phaseDiff<-tol)
        offset += (float)(2 * M_PI);
			else if(phaseDiff>tol)
        offset -= (float)(2 * M_PI);
		}

		uPhaseData[i+1] = phasedata[i+1] + offset;
	}
			
	return;
}

void matrixMultiply(const float* left, int leftRows, int innerDimension, const float* right, int rightColumns, float *output)
{
	int m,n,p;
	for(m=0; m<leftRows; m++)
	{
		for(n=0; n<rightColumns; n++)
		{
			output[m+leftRows*n] = 0;
			for(p=0; p<innerDimension; p++)
				output[m+leftRows*n] += left[m+leftRows*p]*right[p+leftRows*n];
		}
	}
}

float modulus(float b, float m)
{
	return (b/m - (int)(b/m))*m;
}

int leastCommonMultiple(int a, int b)
{
  int l=1;
  for(int n=1; n<=MIN(a,b); n++)
    if((n%a==0) && (n%b==0))
      l = n;

  return l;
}
