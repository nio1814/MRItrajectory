#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <vector>
#include <numeric>

extern "C"
{
#include "cones.h"
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

pybind11::dict loadTrajectory(std::string filepath)
{
  TrajectoryGenerator loader;
  Trajectory* trajectoryLoaded = loader.load(filepath);

  pybind11::dict trajectory;

  pybind11::list fieldOfView;
  pybind11::list spatialResolution;
  for(int d=0; d<trajectoryLoaded->numDimensions; d++)
  {
    fieldOfView.append(trajectoryLoaded->fieldOfView[d]);
    spatialResolution.append(trajectoryLoaded->spatialResolution[d]);
  }
  trajectory["field_of_view"] = fieldOfView;
  trajectory["spatial_resolution"] = spatialResolution;
  trajectory["sampling_interval"] = trajectoryLoaded->samplingInterval;

  std::vector<int> size = {trajectoryLoaded->numReadouts, trajectoryLoaded->numDimensions, trajectoryLoaded->numWaveformPoints};
  trajectory["gradients"] = createNumpyArray(trajectoryLoaded->gradientWaveforms, size);

  size = {trajectoryLoaded->numReadouts, trajectoryLoaded->numDimensions, trajectoryLoaded->numReadoutPoints};
  trajectory["kspace"] =  createNumpyArray(trajectoryLoaded->kSpaceCoordinates, size);

  if(trajectoryLoaded->densityCompensation)
  {
    size = {trajectoryLoaded->numReadouts, trajectoryLoaded->numReadoutPoints};
    trajectory["density_compensation"] =  createNumpyArray(trajectoryLoaded->densityCompensation, size);
  }

  if(trajectoryLoaded->type == CONES)
  {
    Cones* cones = loader.cones();
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
      trajectory["interpolation"] = interpolation;
      deleteCones(&cones);
    }
  }


  /*switch(loader.trajectoryType())
  {
    case CONES:
      break;
    default:
      break;
  }*/


  return trajectory;
}

}

PYBIND11_MODULE(cmri, module)
{
/*  pybind11::class_<CTrajectory>("Trajectory", pybind11::init<std::string>())
      .def_readonly("sampling_interval", &CTrajectory::samplingInterval)
      .def_readonly("num_dimensions", &CTrajectory::numDimensions)
      .def_readonly("num_readouts", &CTrajectory::numReadouts)
      .def_readonly("gradients", &CTrajectory::gradients)
      .def_readonly("kspace", &CTrajectory::kSpace)
      ;*/
  module.def("load_trajectory", &py::loadTrajectory);
}
