#ifndef SPINWARP_H
#define SPINWARP_H

//#include "trajectory_s.h"
//#include "vd.h"

struct PhaseEncoder;
struct Trajectory;
struct VariableDensity;
enum WaveformStorageType;

enum SWorder{swRASTER, swVDRAD};

struct Cartesian
{
  float* readoutGradient;
  int numReadoutPoints;
  int numReadoutWaveformPoints;
  int readoutAxis;
  int numPreReadoutPoints;
  int numPreReadoutWaveformPoints;
  float readoutFieldOfView;
  float readoutSpatialResolution;
  struct PhaseEncoder* phaseEncoders[2];
  float phaseEncodeFieldOfView[2];
  float phaseEncodeSpatialResolution[2];
  int numPrePhaseEncodeWaveformPoints[2];
  int phaseEncodeAxes[2];
  int numPhaseEncodeAxes;
  int numWaveformPoints;
  float samplingInterval;
};

struct CartesianSchedule
{
	int length;
  int *phaseEncodeIndex1;	/**< 1st dimension phase encode scaling */
  int *phaseEncodeIndex2;	/**< 2nd dimension phase encode scaling */
};

struct CartesianSchedule* newCartesianSchedule(const int length, const int numPhaseEncodeAxes);
//void spinwarpInfoToTrajectory(struct SpinwarpInfo *swinfo, struct Trajectory *traj);

//void initSpinwarpInfo(struct SpinwarpInfo *swinfo);
//int spinwarpSaveInfo(const char* filename, struct SpinwarpInfo *info);
//int spinwarpLoadInfo(const char* filename, struct SpinwarpInfo *info);

//struct Cartesian* newCartesian(sched, int length, int ndim);
//int spinwarpSaveSched(const char* filename, struct SpinwarpSched *sched);

struct Cartesian generateCartesian(float* fieldOfView, float *spatialResolution, const char *axisOrder, int numDimensions, int doRephase, float samplingInterval, float gradientLimit, float slewLimit);
struct Trajectory* cartesianToTrajectory(const struct Cartesian* cartesian, const enum WaveformStorageType storage);
//void calcNDFTrdout(int n, struct SpinwarpSched *sched, struct SpinwarpInfo *swinfo, float *gx, float *gy, float *gz);

#endif

