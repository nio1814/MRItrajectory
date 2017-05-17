#ifndef MATHOPS_H
#define MATHOPS_H

int roundToInteger(float num);
void calculatePhase(float *realdata, float *imagdata, float *phasedata, long npts, float phasestart, float phasescale);
void calculateMagnitude(float *realdata, float *imagdata, float *magdata, long npts);
void unwrapPhase(float *phasedata, float* uPhaseData, int numPoints, int frsize, float tol);

#endif // MATHOPS_H
