#include "rings.h"

#include "variabledensity.h"
#include "mathops.h"
#include "arrayops.h"
#include "mrgradient.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef GMR
#define GMR 4257.6388f
#endif

void ccr_calcGT(float gradientLimit, float slewRateLimit, float kR, float fieldOfView, float T, float *G0, float *T0)
{
  float T_halfsine0;

/* Sampling constraints */
float readoutGradientLimit = calculateMaxReadoutGradientAmplitude(fieldOfView, T);
float Gmax = MIN(readoutGradientLimit, gradientLimit);

/* Timing
 % 1/cm, area accrued during half-cycle of sine should be 2*kR
          % with no prewinder, kx would be 0 ~ -2kR in half a cycle
          % with initial offset of kR, kx would be kR ~ -kR*/

float A = 2*kR;

if( fabs(A/GMR)< 2*Gmax*Gmax/slewRateLimit )
    /* slew-limited */
  T_halfsine0 = M_PI*sqrt( fabs(A/GMR)/2/slewRateLimit );
else
    /* grad-limited */
  T_halfsine0 = M_PI*fabs(A/GMR)/2/Gmax;

/* calculate grad amplitude from T_halfsine */
*G0 = M_PI*fabs(A/GMR)/2/T_halfsine0;

/* calculate period */
*T0 = 2*T_halfsine0;

  return;
}


/**
% Project: Concentric Rings
% filename: ccr_gengrads.m
%           generate gradient waveforms
% 03/01/2006
% Hochong Wu

function [ktraj, grad, slew, time, theta_rd, lenStruct] = ...
    ccr_gengrads(SYS, Tw, Ts, FOV, npe, REVLN, USAMP, lenRd_tgt, gdepIN, DBGCCR)
% SYS stucture includes:
% % Constants
% GMR = 4258; % Hz/G
% % Hardware constraints
% TMIN = 4 *10^-6;            % 4 us
% GMAX = 40 *(10/100);        % 40 mT/m --> 4 G/cm
% SMAX = 150 *(10/100/10^-3); % 150 mT/m/ms --> 15,000 G/cm/s
*/
int generateRingGradients(struct VariableDensity* variableDensity, float fieldOfView, float spatialResolution, float readoutDuration, float numRevolutions, float gradientLimit, float slewRateLimit, float samplingInterval, float **gradx, float **grady, int* numDephasePoints, int* numReadoutPoints, int* numRephasePoints, float **theta_rd)
{
  int useLTGT = 1;
  int sf;
  int nsrev2;
  float gradientTimeScaleFactor;
  float *gringx, *gringy;
  float *kringx, *kringy;
  float kxpre;
  float kypre;
  int useGDEPIN = 0;
  int Nrep0;
  float *kdep;
  float *gdep0 = NULL;
  float *sdep;
  float *tdep;

  int n;

  float krMax = 5/spatialResolution;

/*** G READOUT ***
% ------------------------------------------------------------------------
[Gtgt, T0] = ccr_calcGT( GMR,gradMaxSys,slewMaxSys, kR,FOV/USAMP,Ts, DBGCCR );*/
float fieldOfViewInitial;
getFieldOfView(variableDensity, 0, &fieldOfView, &fieldOfViewInitial, 1);

float maxGradientAmplitude;
float minReadoutDuration;
ccr_calcGT(gradientLimit, slewRateLimit, krMax, fieldOfViewInitial, samplingInterval, &maxGradientAmplitude, &minReadoutDuration);
/* number of samples for REVLN*(2pi) = REVLN * period of sine /(sec/samples) */
*numReadoutPoints = ceil(numRevolutions*minReadoutDuration/samplingInterval);
/* make length a multiple of 8 ... for SE sequences
%nsrev = 8*ceil(nsrev/8);
% 2009/10/29 also make length a multiple of 2*REVLN to ensure even#samps/2pi*/
sf = leastCommonMultiple(8, 2*numRevolutions);
*numReadoutPoints = sf*ceil(1.0**numReadoutPoints/sf);
if( useLTGT ) /* align to target readout window duration */
{
  *numReadoutPoints = readoutDuration/samplingInterval;
    nsrev2 = sf*ceil(1.0**numReadoutPoints/sf); /* make sure it's a multiple of (8,2*REVLN) */
    gradientTimeScaleFactor = minReadoutDuration/readoutDuration;
    if( gradientTimeScaleFactor>1 ) /* nsrev2 < nsrev */
  {
        /* this will amplify gradient amplitude! */
        fprintf(stderr, "ccr_gengrads: gradientTimeScaleFactor > 1! Cannot align readout length to specs.\n");
        fprintf(stderr, "ccr_gengrads: Min = %d, Spec = %d\n", *numReadoutPoints, nsrev2);
        fprintf(stderr, "ccr_gengrads: Run script again with lenRdtgt>=%d\n", *numReadoutPoints);
        fprintf(stderr, "ccr_gengrads: END\n");
        return 1;
  }
    else /* nsrev < nsrev2 */
        /* this means scaling down gradient amplitude
        % and stretching the readout window from nsrev to nsrev2*/
        *numReadoutPoints = nsrev2;
}
else
    gradientTimeScaleFactor = 1;

/* % now we know theta for the sinusoids
%thetaRd = REVLN*2*pi*[0:nsrev-1]/nsrev;
%theta_rd = REVLN*2*pi*[0:nsrev-1]/(nsrev-1);

% now we know theta for the sinusoids
% 0 and REVLN*2pi come from the last/first sample of the dephaser/rephaser
%theta_rd = REVLN*2*pi*[1:nsrev]/(nsrev+1);
% 2009/10/29 use an even# of samps over 2pi
theta_rd = REVLN*2*pi*[1:nsrev]/(nsrev);*/
*theta_rd = (float*)malloc(*numReadoutPoints*sizeof(float));
for(n=0; n<*numReadoutPoints; n++)
  (*theta_rd)[n] = n*numRevolutions*2*M_PI/(*numReadoutPoints);

/* use ATTEN to scale the sinusoids
gringx = -ATTEN*Gtgt*sin(theta_rd);
gringy =  ATTEN*Gtgt*cos(theta_rd);
gring  = gringx + i*gringy;*/

gringx = (float*)malloc(*numReadoutPoints*sizeof(float));
gringy = (float*)malloc(*numReadoutPoints*sizeof(float));

for(n=0; n<*numReadoutPoints; n++)
{
  (*theta_rd)[n] = numRevolutions*2*M_PI*(n+1)/(*numReadoutPoints);
  gringx[n] = -gradientTimeScaleFactor*maxGradientAmplitude*sin((*theta_rd)[n]);
  gringy[n] =  gradientTimeScaleFactor*maxGradientAmplitude*cos((*theta_rd)[n]);
}

/* calculate prewinder area
kring = GMR*Tw*cumsum(gring);*/
kringx = (float*)malloc(*numReadoutPoints*sizeof(float));
kringy = (float*)malloc(*numReadoutPoints*sizeof(float));

cumulativeSum(gringx, kringx, *numReadoutPoints);
cumulativeSum(gringy, kringy, *numReadoutPoints);
scalecomplex(kringx, kringy, GMR*samplingInterval, 0.0f, *numReadoutPoints);

/*%figure; plot( real(kring) ); hold on; plot( imag(kring), 'r--' );
% kxpre to place staring point of rings at (kR, 0)*/
/*kxpre = -min(real(kring))/2;*/
kxpre  = -minValue(kringx, *numReadoutPoints)/2;
/*% kypre gets rid of jitters in sine due to discretization and boundary cond'ns
kypre = -(fmaxf(imag(kring))+fminf(imag(kring)))/2;*/
kypre = -(maxValue(kringy, *numReadoutPoints) + minValue(kringy,*numReadoutPoints))/2;


/* *** G DEPHASER/REPHASER ***
% ------------------------------------------------------------------------
% Target area
Ktgt = kxpre + i*kypre;
% Prewinder*/
  int numDephasePointsInitial = 0;
if( useGDEPIN ); /* generate gdep0 based on input gdepIN
     scale gdepx to provide kxpre,
    adepx  = GMR*Tw*sum( real(gdepIN) );
  adepx = GMR*Tw*sumfloats(gdepIN,
    gdepx0 = abs(kxpre/adepx)*real(gdepIN);
    % extract gdepyT and gdepyR
    idx = max( find( imag(gdepIN)==0 ) );
    gdepyT0 = imag(gdepIN(1:idx));
    gdepyR0 = imag(gdepIN(idx:end));
    % scale gdepyR and gdepyT to provide Gtgt,
    gtgt_y  = imag( gdepIN(end) );
    gdepyR0 = abs(ATTEN*Gtgt/gtgt_y)*gdepyR0;
    gdepyT0 = abs(ATTEN*Gtgt/gtgt_y)*gdepyT0;
    % scale gdepyT to provide kypre
    adepyR  = GMR*Tw*sum( gdepyR0 );
    adepyT  = GMR*Tw*sum( gdepyT0 );
    gdepyT0 = abs((-adepyR+kypre)/adepyT)*gdepyT0;
    % now we have gdepy0
    gdepy0 = [gdepyT0 gdepyR0(2:end)];
    % combine gxdep and gydep
    gdep0 = gdepx0 + i*gdepy0;*/
else
  {
    /*[kdep, gdep0, sdep, tdep] = bipolar_prewinder( SYS,Tw, Ktgt,ATTEN*Gtgt, DBGCCR );*/
    float readoutGradientLimit = calculateMaxReadoutGradientAmplitude(fieldOfView, samplingInterval);
  grd_bipolar_prewinder(readoutGradientLimit, gradientLimit, slewRateLimit, samplingInterval, kxpre, kypre, gradientTimeScaleFactor*maxGradientAmplitude, &kdep, &gdep0, &sdep, &tdep, &numDephasePointsInitial);
  }

/* zeropad #samps to a multiple of 4 */
/*gdep = zpad4( gdep0, 1 );*/
  *numDephasePoints = numDephasePointsInitial + numDephasePointsInitial%4;

/* Rewinder
grep0 = -real( fliplr(gdep0) ) + i*imag( fliplr(gdep0) );*/

/* % 2009/10/29 explicit calculation
% KtgtRep = GMR*Tw*sum([gdep gring]);
% [krep, grep0, srep, trep] = bipolar_prewinder( SYS,Tw, KtgtRep,ATTEN*Gtgt, DBGCCR );
% grep0 = -real( fliplr(grep0) ) + i*imag( fliplr(grep0) );

% zeropad #samps to a multiple of 4
% take samples 2:end, since last/first samples of gring/grep are the same
grep = zpad4( grep0(2:end), 0 );*/
Nrep0 = *numDephasePoints-1;
  *numRephasePoints = Nrep0 + Nrep0%4;

/* pad grep and gdep to the same length
if( numel(grep)<numel(gdep) )
    grep = zpad2L( grep, numel(gdep), 0 );
else
    gdep = zpad2L( gdep, numel(grep), 1 );
end */

if(numRephasePoints>numDephasePoints)
  numDephasePoints = numRephasePoints;
else if(numDephasePoints<numRephasePoints)
  numRephasePoints = numDephasePoints;

/* *** POST-PROCESSING ***
% ------------------------------------------------------------------------
% Combine processed gradients to output
grad  = [gdep gring grep];
ktraj = GMR*Tw*cumsum( grad );
slew  = diff([0 grad])/Tw;
time  = [0:numel(grad)-1]*Tw;*/

  int numPoints = *numDephasePoints + *numReadoutPoints + *numRephasePoints;
*gradx = (float*)calloc(numPoints, sizeof(float));
*grady = (float*)calloc(numPoints, sizeof(float));
memcpy(&(*gradx)[*numDephasePoints-numDephasePointsInitial], gdep0, sizeof(float)*numDephasePointsInitial);
memcpy(&(*grady)[*numDephasePoints-numDephasePointsInitial], &gdep0[numDephasePointsInitial], sizeof(float)*numDephasePointsInitial);

memcpy(&(*gradx)[*numDephasePoints], gringx, sizeof(float)**numReadoutPoints);
memcpy(&(*grady)[*numDephasePoints], gringy, sizeof(float)**numReadoutPoints);

for(n=0; n<*numRephasePoints; n++)
{
    (*gradx)[*numDephasePoints+*numReadoutPoints+n] = -(*gradx)[*numDephasePoints-1-n];
    (*grady)[*numDephasePoints+*numReadoutPoints+n] = (*grady)[*numDephasePoints-1-n];
}

/* Rotate rephaser by fractional part of number of revolutions */
scalecomplex(&(*gradx)[*numDephasePoints+*numReadoutPoints], &(*grady)[*numDephasePoints+*numReadoutPoints], cos(2*M_PI*numRevolutions), sin(2*M_PI*numRevolutions), *numRephasePoints);

/* lenStruct
lenStruct.lendep = numel(gdep);  % length of gdep in 4us samples
lenStruct.lenRd  = numel(gring); % number of 4us readout points /ring
lenStruct.lenrep = numel(grep);  % length of grep in 4us samples
lenStruct.lenTot = numel(grad);  % total length of gradient in 4us samples*/

  free(gringx);
  free(gringy);
  free(kringx);
  free(kringy);

  free(kdep);
  free(gdep0);
  free(sdep);
  free(tdep);

        return 0;
}

struct Trajectory *generateRings(struct VariableDensity *variableDensity, float fieldOfView, float spatialResolution, float readoutDuration, float gradientLimit, float slewRateLimit, float samplingInterval)
{
  struct Trajectory* trajectory = newTrajectory();

  /* make matrix size even */
  adjustSpatialResolution(fieldOfView, trajectory->imageDimensions, &spatialResolution);

  trajectory->samplingInterval = samplingInterval;
  for (int n=0; n<2; n++)
  {
    trajectory->fieldOfView[n] = fieldOfView;
    trajectory->spatialResolution[n] = spatialResolution;
  }
//  ccrparams[8] = nrev;
  trajectory->numBases = 1;
  trajectory->type = RINGS;
  trajectory->storage = STORE_ALL;
  trajectory->imageDimensions[1] = trajectory->imageDimensions[0];
  trajectory->numDimensions = 2;

  float* gradientBasis[2];
  float *theta;
  int numRephasePoints;

  generateRingGradients(variableDensity, fieldOfView, spatialResolution, readoutDuration, 1,
                        gradientLimit, slewRateLimit, samplingInterval,
                        &gradientBasis[0], &gradientBasis[1],
                        &trajectory->numPreReadoutPoints, &trajectory->numReadoutPoints, &numRephasePoints, &theta);

//  int numReadouts;
//    if(vd->nsteps==1)
//    {
//      /* full phase encoding */
//      ccr->ninter = ccrparams[7]+1;

//      if(ccr->ninter)
//        free(ccinfo->pes);
//      ccinfo->pes = (float*)calloc(ccr->ninter, sizeof(float));
//      calcPeScales(ccr->ninter, ccinfo->pes, 1, peHALF);
//    }
//    else
//      /* variable step phase encoding */
////    n = ccr->ninter;
//    calcVdPeScales(vd, 0, spatialResolution, 1, peHALF, &ccinfo->pes, &numReadots);

  int numPhaseEncodes = trajectory->imageDimensions[0]/2 + 1;
    allocateTrajectory(trajectory, trajectory->numReadoutPoints, trajectory->numPreReadoutPoints+trajectory->numReadoutPoints+numRephasePoints, 2, 1, numPhaseEncodes, STORE_ALL);

  float* kSpaceWaveform = (float*)malloc(trajectory->numWaveformPoints*sizeof(float));
    for(int r=0; r<trajectory->numReadouts; r++)
  {
      float scale = r/(float)trajectory->numReadouts;
      for(int d=0; d<2; d++)
      {
        float* gradient = trajectoryGradientWaveform(trajectory, r, d);
        memcpy(gradient, gradientBasis[d], trajectory->numWaveformPoints*sizeof(float));
//     scale = ccinfo->pes[n];


        scalefloats(gradient, trajectory->numWaveformPoints, scale);

        gradientToKspace(gradient, kSpaceWaveform, samplingInterval, trajectory->numWaveformPoints);
        memcpy(trajectoryKspaceWaveform(trajectory, r, d), &kSpaceWaveform[trajectory->numPreReadoutPoints], trajectory->numReadoutPoints*sizeof(float));
      }

        float* dcf = trajectoryDensityCompensationWaveform(trajectory, r);
        if(r==0)
            setValues(dcf, trajectory->numReadoutPoints, 1.0f/trajectory->numReadoutPoints);
        else
            setValues(dcf, trajectory->numReadoutPoints, scale);
  }

  return trajectory;
}
