#include "phaseencode.h"

#include "mrgradient.h"
#include "arrayops.h"

#include <stdlib.h>
#include <string.h>

void phaseEncodeScales(float fieldOfView, float spatialResolution, float **scales, float direction, int fullPhaseEncoding, int *numPhaseEncodes)
{
  *numPhaseEncodes = fieldOfView*10/spatialResolution;
  if(!fullPhaseEncoding)
    *numPhaseEncodes /= 2;
  *scales = (float*)malloc(*numPhaseEncodes*sizeof(float));

  int p;

  if(fullPhaseEncoding)
  {
    for(p=0; p<*numPhaseEncodes; p++)
     if(*numPhaseEncodes%2)
      (*scales)[p] = direction*2.0f*(p-*numPhaseEncodes/2)/(float)*numPhaseEncodes;
     else
      (*scales)[p] = direction*(2.0f*p-*numPhaseEncodes)/(*numPhaseEncodes);
  }
  else
  {
      for(p=0; p<*numPhaseEncodes; p++)
        if(direction>0)
      (*scales)[p] = p/(float)*numPhaseEncodes;
      else
          (*scales)[p] = (*numPhaseEncodes-p-1)/(float)*numPhaseEncodes;
  }

  return;
}

struct PhaseEncoder* newPhaseEncoder(float fieldOfView, float spatialResolution, int fullPhaseEncoding, int direction, int rewind, int numPointsBetween, float samplingInterval, float gradientLimit, float slewLimit)
{
  struct PhaseEncoder* phaseEncoder = (struct PhaseEncoder*)malloc(sizeof(struct PhaseEncoder));
  phaseEncoder->fieldOfView = fieldOfView;
  phaseEncoder->spatialResolution = spatialResolution;

  int numPointsRamp;
  float* phaseEncodeGradient;
  grd_trapezoid(gradientLimit, slewLimit, 5/spatialResolution, samplingInterval, NULL, &phaseEncodeGradient, NULL, NULL, &phaseEncoder->numPointsPhaseEncode, &numPointsRamp);

  phaseEncodeScales(fieldOfView, spatialResolution, &phaseEncoder->scales, direction, fullPhaseEncoding, &phaseEncoder->numPhaseEncodes);

  phaseEncoder->numPoints = phaseEncoder->numPointsPhaseEncode + numPointsBetween;
  if(rewind)
    phaseEncoder->numPoints += phaseEncoder->numPointsPhaseEncode;

  phaseEncoder->gradient = (float*)calloc(phaseEncoder->numPoints, sizeof(float));

  memcpy(phaseEncoder->gradient, phaseEncodeGradient, phaseEncoder->numPointsPhaseEncode*sizeof(float));

  if(rewind)
  {
    float* gradientRewind = &phaseEncoder->gradient[phaseEncoder->numPointsPhaseEncode+numPointsBetween];
    memcpy(gradientRewind, phaseEncodeGradient, phaseEncoder->numPointsPhaseEncode*sizeof(float));
    scalefloats(gradientRewind, phaseEncoder->numPointsPhaseEncode, -1);
  }

  free(phaseEncodeGradient);

  return phaseEncoder;
}

void deletePhaseEncoder(struct PhaseEncoder **phaseEncoder)
{
    free((*phaseEncoder)->gradient);
    free(*phaseEncoder);
}
