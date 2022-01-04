#include "trajectorygenerator.h"

extern "C"
{
#include "spiral.h"
#include "radial.h"
#include "cones.h"
#include "rings.h"
#include "spinwarp.h"
}

#include <algorithm>
#include <sstream>


TrajectoryGenerator::TrajectoryGenerator(TrajectoryType type) :
  m_trajectoryType(type),
  m_variableDensity(newVariableDensity())
{
  resetVariableDensity();
}

TrajectoryGenerator::~TrajectoryGenerator()
{
  deleteVariableDensity(&m_variableDensity);
}

TrajectoryType TrajectoryGenerator::trajectoryType()
{
  return m_trajectoryType;
}

void TrajectoryGenerator::setTrajectoryType(TrajectoryType type)
{
  m_trajectoryType = type;
}

void TrajectoryGenerator::setFieldOfView(double* fov, int numDims)
{
  for(int n=0; n<numDims; n++)
    m_fieldOfView[n] = fov[n];
}

void TrajectoryGenerator::setFieldOfView(std::vector<float> fov)
{
  int numDims = std::min(fov.size(), (size_t)3);
  for(int n=0; n<numDims; n++)
    m_fieldOfView[n] = fov[n];
}

void TrajectoryGenerator::setFieldOfView(const float extent, const int axis)
{
  m_fieldOfView[axis] = extent;
}

std::vector<float> TrajectoryGenerator::fieldOfView()
{
  return std::vector<float>(m_fieldOfView, m_fieldOfView + this->numDimensions());
}
float TrajectoryGenerator::maxFieldOfView()
{
  return *std::max_element(m_fieldOfView, m_fieldOfView + numDimensions());
}

float TrajectoryGenerator::maxFieldOfViewXY()
{
  return std::max(m_fieldOfView[0], m_fieldOfView[1]);
}

void TrajectoryGenerator::setFilterFieldOfView(float filterFOV)
{
  m_filterFieldOfView = filterFOV;
}

void TrajectoryGenerator::setSpatialResolution(double *spatialRes, int numDims)
{
  for(int n=0; n<numDims; n++)
    m_spatialResolution[n] = spatialRes[n];

}

void TrajectoryGenerator::setSpatialResolution(std::vector<float> resolution)
{
  std::copy(resolution.data(), resolution.data() + resolution.size(),  m_spatialResolution);
}

void TrajectoryGenerator::resetVariableDensity()
{
  setToFullySampled(m_variableDensity);
}


float TrajectoryGenerator::minSpatialResolution()
{
  return *std::min_element(m_spatialResolution, m_spatialResolution + numDimensions());
}

float TrajectoryGenerator::minSpatialResolutionXY()
{
  return std::min(m_spatialResolution[0], m_spatialResolution[1]);
}

void TrajectoryGenerator::setSamplingInterval(float interval)
{
  m_samplingInterval = interval;
}

void TrajectoryGenerator::setGradientAmplitudeLimit(float gradientLimit)
{
  m_gradientLimit = gradientLimit;
}

void TrajectoryGenerator::setSlewRateLimit(float slewLimit)
{
  m_slewRateLimit = slewLimit;
}


void TrajectoryGenerator::setReadoutDuration(float readoutDuration)
{
  m_readoutDuration = readoutDuration;
}

void TrajectoryGenerator::setNumReadouts(int numReadouts)
{
  m_numReadouts = numReadouts;
}

void TrajectoryGenerator::setRotatable(bool status)
{
  m_rotatable = status;
}

void TrajectoryGenerator::setNumBases(int bases)
{
  m_numBases = bases;
}

void TrajectoryGenerator::setCompensation(const std::string compensation)
{
  m_compensation = compensation;
}

void TrajectoryGenerator::setStorage(WaveformStorageType type)
{
  m_storage = type;
}

WaveformStorageType TrajectoryGenerator::storage()
{
  return m_storage;
}

bool TrajectoryGenerator::generate()
{
  switch(m_trajectoryType)
  {
    case SPIRAL:
      m_trajectory = generateSpirals(m_variableDensity, m_fieldOfView[0], m_spatialResolution[0], m_readoutDuration, 1, m_samplingInterval, m_numReadouts, ARCHIMEDEAN, 0, m_fieldOfView[0], m_gradientLimit, m_slewRateLimit);
      break;
    case RADIAL:
      m_trajectory = generateRadial2D(m_fieldOfView[0], m_fieldOfView[1], EllipticalShape, m_spatialResolution[0], m_spatialResolution[1], EllipticalShape, m_fullProjection, 1, m_gradientLimit, m_slewRateLimit, m_samplingInterval);
      break;
    case RADIAL3D:
      m_trajectory = generateRadial3D(m_fieldOfView[0], m_fieldOfView[1], m_fieldOfView[2],
          EllipticalShape, EllipticalShape,
          m_spatialResolution[0], m_spatialResolution[1], m_spatialResolution[2],
          m_fullProjection, getFinalScale(m_variableDensity), m_gradientLimit, m_slewRateLimit, m_samplingInterval);
      break;
    case CONES:
    {
      if(m_cones)
        deleteCones(&m_cones);

      InterConeCompensation compensation = NO_COMPENSATION;
      if(m_compensation == "none")
        compensation = NO_COMPENSATION;
      else if(m_compensation == "1")
        compensation = Compensation1;
      else if(m_compensation == "2")
        compensation = Compensation2;

      m_cones = generateCones(maxFieldOfViewXY(), m_fieldOfView[2], m_variableDensity, m_spatialResolution[0], m_spatialResolution[2], m_numBases, 1, compensation, m_readoutDuration, m_samplingInterval, maxFieldOfView(), m_gradientLimit, m_slewRateLimit, STORE_ALL);
      if(!m_cones)
        return false;
      m_trajectory = m_cones->trajectory;
    }
      break;
//    case SWIRL:
//      if(m_swirl)
//        freeSwirl(m_swirl);
//      m_swirl = generateSwirls(maxFieldOfView(), maxFieldOfView(), m_variableDensity, minSpatialResolution(), 0, m_readoutDuration, m_samplingInterval, m_gradientLimit, m_slewRateLimit, STORE_ALL);
//      m_trajectory = m_swirl->trajectory;
//      break;
    case RINGS:
      m_trajectory = generateRings(m_variableDensity, maxFieldOfView(), minSpatialResolution(), m_readoutDuration, m_gradientLimit, m_slewRateLimit, m_samplingInterval);
      break;
    case STACK_OF_SPIRALS:
      {
        StackOfSpirals* spirals = generateStackOfSpirals(m_variableDensity, maxFieldOfViewXY(), m_fieldOfView[2], minSpatialResolutionXY(), m_spatialResolution[2], m_readoutDuration, true, m_samplingInterval, 0, m_filterFieldOfView, m_gradientLimit, m_slewRateLimit);
        m_trajectory = stackOfSpiralsToTrajectory(spirals, STORE_ALL);
        deleteStackOfSpirals(&spirals);
      }
      break;
    case CARTESIAN3D:
      const int numDims = m_trajectoryType==CARTESIAN3D ? 3 : 2;
      Cartesian cartesian = generateCartesian(m_fieldOfView, m_spatialResolution, "xyz", numDims, false, m_samplingInterval, m_gradientLimit, m_slewRateLimit);
      m_trajectory = cartesianToTrajectory(&cartesian, STORE_ALL);
      break;
  }

  return m_trajectory;
}

bool TrajectoryGenerator::save(std::string filepath)
{
  bool status = false;
  switch(m_trajectoryType)
  {
    case CONES:
      if(m_cones)
        status = !saveCones(filepath.c_str(), m_cones);
      break;
    default:
      status = !saveTrajectory(filepath.c_str(), m_trajectory);
  }

  return status;
}

Trajectory* TrajectoryGenerator::load(std::string filepath)
{
  Trajectory* trajectory = loadTrajectory(filepath.c_str(), LittleEndian);
  TrajectoryType trajectoryType = trajectory->type;
  switch(trajectoryType)
  {
    case CONES:
      if(m_cones)
        deleteCones(&m_cones);
      m_cones = loadCones(filepath.c_str(), LittleEndian, STORE_ALL);
      deleteTrajectory(&trajectory);
      trajectory = m_cones->trajectory;
      break;
    default:
      break;
  }

  m_trajectory = trajectory;
  setTrajectoryType(trajectory->type);

  return m_trajectory;
}

int TrajectoryGenerator::numDimensions()
{
  int numDims = 0;

  switch(m_trajectoryType)
  {
    case SPIRAL:
    case RADIAL:
    case RINGS:
      numDims = 2;
      break;
    case CONES:
    case RADIAL3D:
    case STACK_OF_SPIRALS:
    case CARTESIAN3D:
      numDims = 3;
      break;
  }

  return numDims;
}

std::string TrajectoryGenerator::info()
{
  std::ostringstream text;

  std::string trajectoryName;
  switch(m_trajectoryType)
  {
    case CONES:
      trajectoryName = "Cones";
      break;
    default:
      break;
  }

  text << trajectoryName << " trajectory" << std::endl;
  text << "Readout duration: " << m_readoutDuration * 1e3 << " ms" << std::endl;

  return text.str();
}

Trajectory *TrajectoryGenerator::trajectory()
{
  return m_trajectory;
}

Cones *TrajectoryGenerator::cones() const
{
  if(m_trajectoryType==CONES)
    return m_cones;
  else
    return nullptr;
}


Cones* TrajectoryGenerator::getConesTrajectory()
{
  return m_cones;
}
