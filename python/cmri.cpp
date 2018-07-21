#include <boost/python.hpp>
#include <boost/python.hpp>
#include <boost/python.hpp>

#include <string>
#include <vector>

using NdArray = boost::python::numeric::array;

extern "C"
{
#include "cones.h"
}
#include "trajectorygenerator.h"

NdArray cToNdArray(float* array, std::vector<int> size)
{
  int numPoints = std::accumulate(size.begin(), size.end(), 1, std::multiplies<int>());
  
  boost::python::list data;
  for (int n=0; n<numPoints; n++)
    data.append(array[n]);
  
  NdArray ndarray(data); 
  
  boost::python::list sizeList;
  for (int dimension : size)
    sizeList.append(dimension);
  ndarray.resize(sizeList);

  return ndarray;
}

class CTrajectory
{
public:

  ~CTrajectory()
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
  
  NdArray gradients()
  {
    std::vector<int> size = {m_trajectory->numReadouts, m_trajectory->numDimensions, m_trajectory->numWaveformPoints};

    return cToNdArray(m_trajectory->gradientWaveforms, size);
  }
  
  NdArray kSpace()
  {
    std::vector<int> size = {m_trajectory->numReadouts, m_trajectory->numDimensions, m_trajectory->numReadoutPoints};

    return cToNdArray(m_trajectory->kSpaceCoordinates, size);
  }

private:
  Trajectory* m_trajectory = NULL;
};

namespace CMRI
{
  boost::python::dict loadTrajectory(std::string filepath)
  {
		TrajectoryGenerator loader;
    Trajectory* trajectoryLoaded = loader.load(filepath);

		boost::python::dict trajectory;

		boost::python::list fieldOfView;
		boost::python::list spatialResolution;
		for(int d=0; d<trajectoryLoaded->numDimensions; d++)
		{
			fieldOfView.append(trajectoryLoaded->fieldOfView[d]);
			spatialResolution.append(trajectoryLoaded->spatialResolution[d]);
		}
		trajectory["field_of_view"] = fieldOfView;
		trajectory["spatial_resolution"] = spatialResolution;
    trajectory["sampling_interval"] = trajectoryLoaded->samplingInterval;

    std::vector<int> size = {trajectoryLoaded->numReadouts, trajectoryLoaded->numDimensions, trajectoryLoaded->numWaveformPoints};
    trajectory["gradients"] = cToNdArray(trajectoryLoaded->gradientWaveforms, size);

    size = {trajectoryLoaded->numReadouts, trajectoryLoaded->numDimensions, trajectoryLoaded->numReadoutPoints};
    trajectory["kspace"] =  cToNdArray(trajectoryLoaded->kSpaceCoordinates, size);

		if(trajectoryLoaded->densityCompensation)
		{
			size = {trajectoryLoaded->numReadouts, trajectoryLoaded->numReadoutPoints};
			trajectory["density_compensation"] =  cToNdArray(trajectoryLoaded->densityCompensation, size);
		}
		
    if(trajectoryLoaded->type == CONES)
    {
      Cones* cones = loader.cones();
      if(cones && cones->interpolation)
      {
        boost::python::list readout;
        boost::python::list basis;
        boost::python::list coneIndex;
        boost::python::list scaleXY;
        boost::python::list scaleZ;
        boost::python::list theta;
        boost::python::list thetaIndex;
        boost::python::list phi;
        boost::python::list numInterleavesOnCone;
        boost::python::list interleafOnCone;
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
        boost::python::dict interpolation;
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

BOOST_PYTHON_MODULE(cmri)
{
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

/*  boost::python::class_<CTrajectory>("Trajectory", boost::python::init<std::string>())
      .def_readonly("sampling_interval", &CTrajectory::samplingInterval)
      .def_readonly("num_dimensions", &CTrajectory::numDimensions)
      .def_readonly("num_readouts", &CTrajectory::numReadouts)
      .def_readonly("gradients", &CTrajectory::gradients)
      .def_readonly("kspace", &CTrajectory::kSpace)
      ;*/
  boost::python::def("load_trajectory", &CMRI::loadTrajectory);
}
