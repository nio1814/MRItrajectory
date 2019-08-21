/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#ifndef TRAJECTORY_H
#define TRAJECTORY_H

/**
 \defgroup	Trajectory	k-Space Trajectory
 @{
*/

#include "variabledensity.h"
#include "convertendian.h"

#include <stdio.h>

enum WaveformStorageType {STORE_BASIS, STORE_ALL};
enum TrajectoryType {SPIRAL, RADIAL, RADIAL3D, CONES, RINGS};

struct Trajectory
{
  int numDimensions;
	int imageDimensions[3];
	float spatialResolution[3];
	float fieldOfView[3];
  int numReadouts;
  int numBases;
	float maxGradientAmplitude;
	float maxReadoutGradientAmplitude;
	float maxSlewRate;
	int numWaveformPoints;
  int numPreReadoutPoints;
  int numReadoutPoints;
	float samplingInterval;
	float *densityCompensation;
	short *gradientWaveformsShort;
	float *gradientWaveforms;
	float *kSpaceCoordinates;
	enum WaveformStorageType storage;
	struct VariableDensity *variableDensity;
  enum TrajectoryType type;
};

struct Trajectory* newTrajectory();

/*!
 * \brief Delete a trajectory from memory
 * \param trajectory  Trajectory to delete
 */
void deleteTrajectory(struct Trajectory** trajectory);

void adjustSpatialResolution(float fieldOfView, int *imageDimension, float *spatialResolution);

void gradientToKspace(float* gradient, float* k, float samplingInterval, int length);

void allocateTrajectory(struct Trajectory *trajectory, int numReadoutPoints, int numWaveformPoints, int dimensions, int numBases, int numReadouts, enum WaveformStorageType storage);

void traverseKspaceFromWaveform(float *gradientOriginalX, float *gradientOriginalY, float *gradientOriginalZ, int pointsOriginal, float *kSpaceCoordinatesFinal, float samplingInterval, float maxGradientAmplitude, float maxSlewRate, float **gradientRewoundX, float**gradientRewoundY, float **gradientRewoundZ, int *pointsRewound);

void traverseKspaceToZero(float *gradientOriginalX, float *gradientOriginalY, float *gradientOriginalZ, int pointsOriginal, float samplingInterval, float maxGradientAmplitude, float maxSlewRate, float **gradientRewoundX, float**gradientRewoundY, float **gradientRewoundZ, int *pointsRewound);

void gradientWaveformToShort(const float *gradientsFloat, int points, float maxGradientAmplitudeScanner, short* gradientsShort);

int writeTrajectory(FILE* file, const struct Trajectory* trajectory);
int saveTrajectory(const char* filename, const struct Trajectory* trajectory);

struct Trajectory* readTrajectory(FILE* file, enum Endian endian);
struct Trajectory* loadTrajectory(const char *filename, enum Endian endian);
enum TrajectoryType loadTrajectoryType(const char *filename);

int saveGradientWaveforms(const char *filename, const float* grad, short dimensions, short interleaves, short points, int numReadoutPoints, float FOV, float maxGradientAmplitude, float maxGradientAmplitudeScanner, float samplingInterval, const char* description, enum Endian endian);

void trajectoryCoordinates(int readoutPoint, int readout, const struct Trajectory *trajectory, float *coordinates, float *densityCompensation);

void setTrajectoryPoint(int readoutPoint, int readout, struct Trajectory *trajectory, const float *coordinates, float densityCompensation);

float* trajectoryGradientWaveform(const struct Trajectory* trajectory, int readout, int axis);

float *trajectoryKspaceWaveform(const struct Trajectory *trajectory, int readout, int axis);

/*!
 * \brief Get the density compensation for a single readout
 * \param trajectory
 * \param readout Readout index
 */
float *trajectoryDensityCompensationWaveform(const struct Trajectory* trajectory, int readout);
float calculateMaxReadoutGradientAmplitude(float fieldOfView, float samplingInterval);

void rotateBasis(float* gxBasis, float* gyBasis, struct Trajectory* trajectory, float angleRange);

#endif

