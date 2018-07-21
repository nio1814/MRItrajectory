#ifndef RINGS_H
#define RINGS_H

#include "trajectory.h"

struct Trajectory* generateRings(struct VariableDensity *variableDensity, float fieldOfView, float spatialResolution, float readoutDuration, float gmaxsys, float slewRateLimit, float samplingInterval);

#endif // RINGS_H
