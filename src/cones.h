#ifndef CONES_V4_H
#define CONES_V4_H

#include "trajectory.h"
#include "variabledensity.h"

struct ConesInterpolation
{
	int readouts;
	int *readout;
	int *basis;
	int *cone;
	float *scaleXY;
	float *scaleZ;
	float *theta;
	int *thetaIndex;
	float *phi;
	int *interleavesOnCone;
	int *interleafOnCone;
};

enum InterConeCompensation{NoCompensation, Compensation1, Compensation2};

struct Cones
{
	struct Trajectory trajectory;
	int rotatable;
	enum InterConeCompensation interconeCompensation;
	float* coneAngles;
	float* coneAngleDensityCompensation;
	int numCones;
	float* basisConeAngles;
	int *basisReadoutPoints;
	int *basisWaveformPoints;
	float *basisGradientWaveforms;
	float *basisKspaceCoordinates;
	struct ConesInterpolation interpolation;
};

struct Cones* allocateCones(int bases);
void freeCones(struct Cones *cones);

struct Cones* generateCones(float xyFieldOfView, float zFieldOfView, const struct VariableDensity *variableDensity, float xySpatialResolution, float zSpatialResolution, int bases, int rotatable, enum InterConeCompensation interConeCompensation, float readoutDuration, float samplingInterval, float filterFieldOfView, float maxGradientAmplitude, float maxSlewRate, enum WaveformStorageType storage);
float* conesBasisGradientWaveform(const struct Cones* cones, int basis, int axis);

#endif

