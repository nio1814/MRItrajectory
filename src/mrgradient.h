#ifndef MRGRADIENT_H
#define MRGRADIENT_H

void readoutGradient(int doDephase, int doRephase, int nread, float res, float gradientLimit, float smax, float Ts, float **g, int *nramp, int *ndep, int *npts);

void spoilerGradient(float gradientLimit, float slewRateLimit, float gradientInitial, float deltaKspace, float gradientFinal, float samplingInterval, float** gradientWaveform, int* points);

#endif // MRGRADIENT_H
