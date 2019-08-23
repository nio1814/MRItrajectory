#ifndef PHASEENCODE_H
#define PHASEENCODE_H

struct PhaseEncoder
{
  float* gradient;
  float* scales;
  int numPhaseEncodes;
  int numPoints;
  int numPointsPhaseEncode;
  float fieldOfView;
  float spatialResolution;
};

struct PhaseEncoder* newPhaseEncoder(float fieldOfView, float spatialResolution, int fullPhaseEncoding, int direction, int rewind, int numPointsBetween, float samplingInterval, float gradientLimit, float slewLimit);

#endif // PHASEENCODE_H
