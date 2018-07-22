#ifndef MRGRADIENT_H
#define MRGRADIENT_H

void grd_readout(int doDephase, int doRephase, int nread, float res, float gmax, float smax, float Ts, float **g, int *nramp, int *ndep, int *npts);

void grd_bipolar_prewinder(float gradMax, float gradMaxSys, float slewMax, float T, float kxpre, float kypre, float Gtgt, float **ktraj, float **grad, float **slew, float **time, int *N);

#endif // MRGRADIENT_H
