#ifndef CONES_V4_H
#define CONES_V4_H

#include "trajectory.h"
#include "variabledensity.h"

struct ConesInterpolation
{
  int numReadouts;
	int *readout;
	int *basis;
	int *cone;
	float *scaleXY;
	float *scaleZ;
	float *theta;
	int *thetaIndex;
	float *phi;
  int *numInterleavesOnCone;
	int *interleafOnCone;
};

enum InterConeCompensation{NO_COMPENSATION, Compensation1, Compensation2};

struct Cones
{
  struct Trajectory* trajectory;
	int rotatable;
	enum InterConeCompensation interconeCompensation;
	float* coneAngles;
	float* coneAngleDensityCompensation;
  int numBases;
  int numCones;
	float* basisConeAngles;
	int *numBasisReadoutPoints;
	int *numBasisWaveformPoints;
  int *numBasisInterleaves;
	float *basisGradientWaveforms;
	float *basisKspaceCoordinates;
  struct ConesInterpolation* interpolation;
};

int saveCones(const char* filename, const struct Cones* cones);
struct Cones* loadCones(const char* filename, enum Endian endian, enum WaveformStorageType storage);

struct Cones* newCones(const int numBases);
void deleteCones(struct Cones **cones);

struct Cones* generateCones(float xyFieldOfView, float zFieldOfView, const struct VariableDensity *variableDensity, float xySpatialResolution, float zSpatialResolution, int bases, int rotatable, enum InterConeCompensation interConeCompensation, float readoutDuration, float samplingInterval, float filterFieldOfView, float maxGradientAmplitude, float maxSlewRate, enum WaveformStorageType storage);
float* conesBasisGradientWaveform(const struct Cones* cones, int basis, int axis);

#endif

