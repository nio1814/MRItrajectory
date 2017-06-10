#ifndef DENSITYCOMPENSATION2D_H
#define DENSITYCOMPENSATION2D_H

struct Trajectory;

void convolutionDensityCompensation(const struct Trajectory *trajectory, int padfront, int* Nall, int numLobes, int numIter, float *W);

#endif // DENSITYCOMPENSATION2D_H
