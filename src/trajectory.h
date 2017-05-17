#ifndef TRAJECTORY_H
#define TRAJECTORY_H

/**
 \defgroup	Trajectory	k-Space Trajectory
 @{
*/

#include "variabledensity.h"
#include "convertendian.h"

enum WaveformStorageType {StoreBasis, StoreAll};

struct Trajectory
{
	int dimensions;
	int imageDimensions[3];
	float spatialResolution[3];
	float fieldOfView[3];
	int readouts;
	int bases;
	float maxGradientAmplitude;
	float maxReadoutGradientAmplitude;
	float maxSlewRate;
	int waveformPoints;
	int readoutPoints;
	float samplingInterval;
	float *densityCompensation;
	short *gradientWaveformsShort;
	float *gradientWaveforms;
	float *kSpaceCoordinates;
	enum WaveformStorageType storage;
	struct VariableDensity *variableDensity;
};

void initializeTrajectory(struct Trajectory* trajectory);

void adjustSpatialResolution(float fieldOfView, int *imageDimension, float *spatialResolution);

void gradientToKspace(float* gradient, float* k, float samplingInterval, int length);

void allocateTrajectory(struct Trajectory *trajectory, int readoutPoints, int waveformPoints, int dimensions, int bases, int readouts, enum WaveformStorageType storage);

void traverseKspaceFromWaveform(float *gradientOriginalX, float *gradientOriginalY, float *gradientOriginalZ, int pointsOriginal, float *kSpaceCoordinatesFinal, float samplingInterval, float maxGradientAmplitude, float maxSlewRate, float **gradientRewoundX, float**gradientRewoundY, float **gradientRewoundZ, int *pointsRewound);

void traverseKspaceToZero(float *gradientOriginalX, float *gradientOriginalY, float *gradientOriginalZ, int pointsOriginal, float samplingInterval, float maxGradientAmplitude, float maxSlewRate, float **gradientRewoundX, float**gradientRewoundY, float **gradientRewoundZ, int *pointsRewound);

void gradientWaveformToShort(const float *gradientsFloat, int points, float maxGradientAmplitudeScanner, short* gradientsShort);

int saveGradientWaveforms(const char *filename, const float* grad, short dimensions, short interleaves, short points, int readoutPoints, float FOV, float maxGradientAmplitude, float maxGradientAmplitudeScanner, float samplingInterval, const char* description, enum Endian endian);

#endif

