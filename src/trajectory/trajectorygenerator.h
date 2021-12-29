#ifndef TRAJECTORYGENERATOR_H
#define TRAJECTORYGENERATOR_H

#include <string>
#include <vector>

extern "C"
{
#include "trajectory.h"
}

struct Cones;
struct Swirl;
struct StackOfSpirals;

class TrajectoryGenerator
{
public:
  TrajectoryGenerator(TrajectoryType type=SPIRAL);
  ~TrajectoryGenerator();

  TrajectoryType trajectoryType();

  void setTrajectoryType(TrajectoryType type);
  void setFieldOfView(double* fov, int numDims);
  void setFieldOfView(std::vector<float> fov);
  std::vector<float> fieldOfView();
  float maxFieldOfView();
  float maxFieldOfViewXY();
  void setFilterFieldOfView(float filterFOV);

  /*!
   * \brief Reset the variable density sampling to fully sampled
   */
  void resetVariableDensity();
  void setSpatialResolution(double* spatialRes, int numDims);
  void setSpatialResolution(std::vector<float> resolution);
  float minSpatialResolution();
  float minSpatialResolutionXY();
  void setSamplingInterval(float interval);
  void setGradientAmplitudeLimit(float gradientLimit);
  void setSlewRateLimit(float slewLimit);
  void setReadoutDuration(float readoutDuration);
  void setNumReadouts(int numReadouts);
  void setRotatable(bool status);
  void setNumBases(int bases);
  void setStorage(WaveformStorageType type);
  WaveformStorageType storage();

  bool generate();
  bool save(std::string filepath);
  Trajectory *load(std::string filepath);
  Trajectory* trajectory();
  Cones* cones() const;
  Cones* getConesTrajectory();

  int numDimensions();

  std::string info();
protected:
  TrajectoryType m_trajectoryType;
  float m_fieldOfView[3] = {28,28,28};
  float m_filterFieldOfView = 28;
  float m_spatialResolution[3] = {2,2,2};
  float m_samplingInterval = 4e-6;
  float m_readoutDuration = 8e-3;
  float m_gradientLimit = 4;
  float m_slewRateLimit = 15000;
  int m_numReadouts;
  bool m_rotatable = true;
  int m_numBases = 32;
  bool m_fullProjection = true;
  WaveformStorageType m_storage = STORE_ALL;

  Trajectory* m_trajectory = 0;
  VariableDensity* m_variableDensity = nullptr;
  Cones* m_cones = 0;
  Swirl* m_swirl = 0;
  StackOfSpirals* m_stackOfSpirals = nullptr;
};

#endif // TRAJECTORYGENERATOR_H
