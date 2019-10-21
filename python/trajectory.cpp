#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <vector>
#include <numeric>

extern "C"
{
#include "cones.h"
#include "spiral.h"
}
#include "trajectorygenerator.h"


namespace py {

pybind11::array_t<float> createNumpyArray(float* array, std::vector<int> size)
{
  pybind11::capsule capsule = pybind11::capsule(array, [](void* data)
    {
      free(data);
    }
  );

  return pybind11::array_t<float>(size, array, capsule);
}


class PyTrajectory
{
public:

  ~PyTrajectory()
  {
    deleteTrajectory(&m_trajectory);
  }

	float samplingInterval()
	{
    return m_trajectory->samplingInterval;
	}
	
	int numDimensions()
	{
    return m_trajectory->numDimensions;
	}

  int numReadouts()
  {
    return m_trajectory->numReadouts;
  }
  
  pybind11::array_t<float> gradients()
  {
    std::vector<int> size = {m_trajectory->numReadouts, m_trajectory->numDimensions, m_trajectory->numWaveformPoints};

    return createNumpyArray(m_trajectory->gradientWaveforms, size);
  }
  
  pybind11::array_t<float> kSpace()
  {
    std::vector<int> size = {m_trajectory->numReadouts, m_trajectory->numDimensions, m_trajectory->numReadoutPoints};

    return createNumpyArray(m_trajectory->kSpaceCoordinates, size);
  }

private:
  Trajectory* m_trajectory;
};

pybind11::dict trajectoryToDict(const Trajectory* trajectory, const std::unique_ptr<TrajectoryGenerator> loader=nullptr)
{
  pybind11::dict pyTrajectory;

  pybind11::list fieldOfView;
  pybind11::list spatialResolution;
  for(int d=0; d<trajectory->numDimensions; d++)
  {
    fieldOfView.append(trajectory->fieldOfView[d]);
    spatialResolution.append(trajectory->spatialResolution[d]);
  }
  pyTrajectory["field_of_view"] = fieldOfView;
  pyTrajectory["spatial_resolution"] = spatialResolution;
  pyTrajectory["sampling_interval"] = trajectory->samplingInterval;

  std::vector<int> size = {trajectory->numReadouts, trajectory->numDimensions, trajectory->numWaveformPoints};
  pyTrajectory["gradients"] = createNumpyArray(trajectory->gradientWaveforms, size);

  size = {trajectory->numReadouts, trajectory->numDimensions, trajectory->numReadoutPoints};
  pyTrajectory["kspace"] =  createNumpyArray(trajectory->kSpaceCoordinates, size);

  if(trajectory->densityCompensation)
  {
    size = {trajectory->numReadouts, trajectory->numReadoutPoints};
    pyTrajectory["density_compensation"] =  createNumpyArray(trajectory->densityCompensation, size);
  }

  if(trajectory->type == CONES)
  {
    Cones* cones = loader->cones();
    if(cones && cones->interpolation)
    {
      pybind11::list readout;
      pybind11::list basis;
      pybind11::list coneIndex;
      pybind11::list scaleXY;
      pybind11::list scaleZ;
      pybind11::list theta;
      pybind11::list thetaIndex;
      pybind11::list phi;
      pybind11::list numInterleavesOnCone;
      pybind11::list interleafOnCone;
      for(int n=0; n<cones->interpolation->numReadouts; n++)
      {
        readout.append(cones->interpolation->readout[n]);
        basis.append(cones->interpolation->basis[n]);
        coneIndex.append(cones->interpolation->cone[n]);
        scaleXY.append(cones->interpolation->scaleXY[n]);
        scaleZ.append(cones->interpolation->scaleZ[n]);
        theta.append(cones->interpolation->theta[n]);
        thetaIndex.append(cones->interpolation->thetaIndex[n]);
        phi.append(cones->interpolation->phi[n]);
        numInterleavesOnCone.append(cones->interpolation->numInterleavesOnCone[n]);
        interleafOnCone.append(cones->interpolation->interleafOnCone[n]);
      }
      pybind11::dict interpolation;
      interpolation["readout"] = readout;
      interpolation["basis"] = basis;
      interpolation["cone_index"] = coneIndex;
      interpolation["scale_xy"] = scaleXY;
      interpolation["scale_z"] = scaleZ;
      interpolation["theta"] = theta;
      interpolation["theta_index"] = thetaIndex;
      interpolation["phi"] = phi;
      interpolation["num_interleaves_on_cone"] = numInterleavesOnCone;
      interpolation["interleaf_on_cone"] = interleafOnCone;
      pyTrajectory["interpolation"] = interpolation;
      deleteCones(&cones);
    }
  }

  return pyTrajectory;
}

pybind11::dict loadTrajectory(std::string filepath)
{
  std::unique_ptr<TrajectoryGenerator> loader = std::make_unique<TrajectoryGenerator>();
  Trajectory* trajectory = loader->load(filepath);

  return trajectoryToDict(trajectory);
}

pybind11::dict pyGenerateSpirals(float fieldOfView, float spatialResolution, float readoutDuration, bool rewindTrajectory, float samplingInterval, int numInterleaves, float readoutFieldOfView, float gradientLimit, float slewLimit)
{
  const SpiralType spiralType = ARCHIMEDEAN;
  const float floretAngle = 0;
  if(!readoutFieldOfView)
    readoutFieldOfView = fieldOfView;
  const Trajectory* trajectory = generateSpirals(NULL, fieldOfView, spatialResolution, readoutDuration, rewindTrajectory, samplingInterval, numInterleaves, spiralType, floretAngle, readoutFieldOfView, gradientLimit, slewLimit);

  return trajectoryToDict(trajectory);
}

pybind11::dict pyGenerateCones(float fieldOfViewXY, float fieldOfViewZ, float spatialResolutionXY, float spatialResolutionZ, float readoutDuration, int numBases, bool rotatable, float samplingInterval, float filterFieldOfView, float maxGradientAmplitude, float maxSlewRate)
{
  Cones* cones = generateCones(fieldOfViewXY, fieldOfViewZ, nullptr, spatialResolutionXY, spatialResolutionZ, numBases, rotatable, InterConeCompensation::NO_COMPENSATION, readoutDuration, samplingInterval, filterFieldOfView, maxGradientAmplitude, maxSlewRate, STORE_ALL);

  return trajectoryToDict(cones->trajectory);
}

}

PYBIND11_MODULE(_trajectory, module)
{
/*  pybind11::class_<CTrajectory>("Trajectory", pybind11::init<std::string>())
      .def_readonly("sampling_interval", &CTrajectory::samplingInterval)
      .def_readonly("num_dimensions", &CTrajectory::numDimensions)
      .def_readonly("num_readouts", &CTrajectory::numReadouts)
      .def_readonly("gradients", &CTrajectory::gradients)
      .def_readonly("kspace", &CTrajectory::kSpace)
      ;*/
  module.def("load_trajectory", &py::loadTrajectory);

  const float gradientLimitDefault = 4;
  const float slewLimitDefault = 15e3;
  const float samplingIntervalDefault = 4e-6;
  module.def("generate_spiral", &py::pyGenerateSpirals,
             pybind11::arg("fieldOfView"), pybind11::arg("spatialResolution"),
             pybind11::arg("readoutDuration"),
             pybind11::arg("rewindTrajectory") = true,
             pybind11::arg("samplingInterval") = samplingIntervalDefault,
             pybind11::arg("numInterleaves") = 0,
             pybind11::arg("readoutFieldOfView") = 0,
             pybind11::arg("gradientLimit") = gradientLimitDefault,
             pybind11::arg("slewLimit") = slewLimitDefault);
  module.def("generate_cones", &py::pyGenerateCones,
             pybind11::arg("fieldOfViewXY"), pybind11::arg("fieldOfViewZ"),
             pybind11::arg("spatialResolutionXY"), pybind11::arg("spatialResolutionZ"),
             pybind11::arg("readoutDuration"),
             pybind11::arg("numBases") = 0,
             pybind11::arg("rotatable") = true,
             pybind11::arg("samplingInterval") = samplingIntervalDefault,
             pybind11::arg("filterFieldOfView") = 0,
             pybind11::arg("gradientLimit") = gradientLimitDefault,
             pybind11::arg("slewLimit") = slewLimitDefault);
}
