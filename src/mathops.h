/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#ifndef MATHOPS_H
#define MATHOPS_H

int roundToInteger(float num);
void calculatePhase(float *realdata, float *imagdata, float *phasedata, long npts, float phasestart, float phasescale);
void calculateMagnitude(float *realdata, float *imagdata, float *magdata, long npts);
void unwrapPhase(float *phasedata, float* uPhaseData, int numPoints, int frsize, float tol);

#endif // MATHOPS_H
