#include "mex.h"

#include "trajectorygenerator.h"

extern "C"
{
#include "cones.h"
//#include "mathfunctions.h"
//#include "bh_array.h"
//#include "vd.h"
}

#include <string>
#include <vector>
#include <map>

#include <string.h>

void copyToStruct(mxArray* mxStruct, const char* fieldName, void* data, int length, mxClassID dataType, mxComplexity complexity)
{
  mxAddField(mxStruct, fieldName);
  mxSetField(mxStruct, 0, fieldName, mxCreateNumericMatrix(length, 1, dataType, complexity));
  void* structData = mxGetData(mxGetField(mxStruct, 0, fieldName));
  if(!structData)	
  { 
    char message[128];
    sprintf(message, "Failed to get field %s\n", fieldName);
    mexErrMsgTxt(message);
  }

  int dataPointSize = 0;
  switch(dataType)
  {
    case mxSINGLE_CLASS:
      dataPointSize = sizeof(float);
    case mxINT32_CLASS:
      dataPointSize = sizeof(int);
      break;
    default:
      break;
  }
  memcpy(structData, data, length*dataPointSize);
}

void mexFunction(int numOutputs, mxArray *outputs[],
  int numInputs, const mxArray *inputs[])
{
  mxAssert(numOutputs, "No outputs specified");
	double *dptr;	/* double pointer for input data */
	mxArray *structField;
	int numParams;

  const char trajectoryFieldName[] = "kSpaceCoordinates";
  const char densityCompensationFieldName[] = "dcf";
    const char gradientsFieldName[] = "gradients";
  const char imageDimensionsFieldName[] = "imageDimensions";
  std::vector<const char*> trajectoryStructFields= {
    imageDimensionsFieldName,
    gradientsFieldName,
    trajectoryFieldName ,
    densityCompensationFieldName
  };

  std::vector<const char*> interpolationStructFields = {
    "index",
    "scaleXY",
    "scaleZ"
  };
	
	
	char msg[128];
  VariableDensity* variableDensity = newVariableDensity();
  addLinearVariableDensityStep(variableDensity, 0, 1);

  mxAssert(numInputs>0, "No trajectory type not given");
  mxAssert(mxIsStruct(inputs[1]), "Parameters not given as struct");
	
  int stringLength = mxGetNumberOfElements(inputs[0])+1;
	std::vector<char> charInput(stringLength);
  mxGetString(inputs[0], charInput.data(), stringLength);
	std::string trajectoryTypeString = std::string(charInput.data());
	std::map<std::string, TrajectoryType> stringToTrajectoryType = {
		{"cones", CONES},
	};

  TrajectoryType trajectoryType = stringToTrajectoryType[trajectoryTypeString];
  TrajectoryGenerator generator(trajectoryType);

  switch(generator.trajectoryType())
  {
    case CONES:
      trajectoryStructFields.push_back("basisReadoutLengths");
    break;
  default:
    break;
  }

  bool fieldOfViewSet = false;
  bool spatialResolutionSet = false;

	if(!stringToTrajectoryType.count(trajectoryTypeString))
	{
		sprintf(msg, "Invalid trajectory type %s\n", trajectoryTypeString.c_str());
		mexErrMsgTxt(msg);
	}
    numParams = mxGetNumberOfFields(inputs[1]);
    for(int n=0; n<numParams; n++)
		{
      structField = mxGetFieldByNumber(inputs[1], 0, n);
      std::string structFieldName = std::string(mxGetFieldNameByNumber(inputs[1],n));
			dptr = mxGetPr(structField);
			
			/*Initial FOV*/
      if(structFieldName == "fieldOfView")
			{
        generator.setFieldOfView(dptr, mxGetNumberOfElements(structField));
        fieldOfViewSet = true;
			}
			/*Final FOV*/
//      else if(structFieldName == "finalFieldOfView")
//				isVD = 1;
//				for(d=0; d<mxGetNumberOfElements(structField); d++)
//					traj.FOV[d+3] = dptr[d];
			/* single or vd fov */
      else if(structFieldName == "fov")
			{
//				if(strcmp(mxGetClassName(structField), "VDfov"))
//				{
//					traj.naxes =mxGetNumberOfElements(structField);
//					for(d=0; d<traj.naxes; d++)
//						traj.FOV[d] = dptr[d];
//				}
//				else
//				{
//					fovClassInput = 1;
//					fovClassInputIdx = n;
//					vdField = mxGetProperty(structField, 0, "fov");
//					traj.naxes =mxGetNumberOfElements(vdField);
//					dptr = mxGetPr(vdField);
//					for(d=0; d<traj.naxes; d++)
//						traj.FOV[d] = dptr[d];
//				}
			}
      else if(structFieldName == "spatialResolution")
      {
        generator.setSpatialResolution(dptr, mxGetNumberOfElements(structField));
        spatialResolutionSet = true;
      }
      /* Sampling interval */
      else if(structFieldName == "dt" || structFieldName == "samplingInterval")
          generator.setSamplingInterval(*dptr);
			/* Maximum gradient amplitude */
      else if(structFieldName == "gmax" || structFieldName=="gradientAmplitudeLimit")
          generator.setGradientAmplitudeLimit(mxGetScalar(structField));
			/* Maximum slew rate */
      else if(structFieldName == "smax" || structFieldName=="slewRateLimit")
          generator.setSlewRateLimit(mxGetScalar(structField));
      else if(structFieldName=="readoutDuration")
        generator.setReadoutDuration(mxGetScalar(structField));
			/*Total number of trajectory points per interleaf*/
//      else if(structFieldName == "npts")
//					traj.npts = *dptr;
//      else if(structFieldName == "vdStart"))
//					traj.vdStart = *dptr;
//      else if(structFieldName == "vdEnd"))
//					traj.vdEnd = *dptr;
//      else if(structFieldName == "vdType"))
//					vdType = (enum VdFunc)rnd(*dptr);
//      else if(structFieldName == "vdParam"))
//					vdParam = *dptr;
//      else if(structFieldName == "readOrder"))
//					readOrder = *dptr;
      else if(structFieldName == "nreadouts"
              || structFieldName=="numReadouts")
//					traj.ninter = *dptr;
        generator.setNumReadouts(mxGetScalar(structField));
//      else if(structFieldName == "nintl"))
//					traj.ninter = mxGetScalar(structField);
//      else if(structFieldName == "ic"))
//					doIC = (int)mxGetScalar(structField);
      else if(structFieldName == "rotatable")
          generator.setRotatable(mxGetScalar(structField));
      else if(structFieldName == "nbasis" || structFieldName == "numBases")
          generator.setNumBases(mxGetScalar(structField));
      else if(structFieldName == "fovfilt" || structFieldName=="filterFieldOfView")
        generator.setFilterFieldOfView(mxGetScalar(structField));
			else
			{
        sprintf(msg, "Unrecognized parameter %s\n", structFieldName.c_str());
				mexErrMsgTxt(msg);
			}
		}

  if(!fieldOfViewSet)
    mexErrMsgTxt("Field of view not set");
	
  if(!spatialResolutionSet)
    mexErrMsgTxt("Spatial resolution not set");


//	if(traj.vdStart>=0 && traj.vdEnd<0)
//		traj.vdEnd = traj.kmax;
	
	/* Setup FOV struct */
//	if(fovClassInput)
//	{
//		structField = mxGetFieldByNumber(inputs[1], 0, fovClassInputIdx);
		
//		initVd(&vd, traj.FOV, traj.naxes, 0);
//		vdField = mxGetProperty(structField, 0, "nsteps");
//		vd.nsteps = mxGetScalar(vdField);
		
//		vdField = mxGetProperty(structField, 0, "kr");
//		dptr = mxGetPr(vdField);
//		for(n=0; n<mxGetNumberOfElements(vdField); n++)
//			vd.kr[n] = dptr[n];
			
//		vdField = mxGetProperty(structField, 0, "scale");
//		dptr = mxGetPr(vdField);
//		for(n=0; n<mxGetNumberOfElements(vdField); n++)
//			vd.scale[n] = dptr[n];
			
//		vdField = mxGetProperty(structField, 0, "param");
//		dptr = mxGetPr(vdField);
//		for(n=0; n<mxGetNumberOfElements(vdField); n++)
//			vd.param[n] = dptr[n];
			
//		vdField = mxGetProperty(structField, 0, "type");
//		dptr = mxGetPr(vdField);
//		for(n=0; n<mxGetNumberOfElements(vdField); n++)
//			vd.type[n] = (enum VdFunc)dptr[n];
//	}
//	else
//	{
//		/* can currently only handle isoptric undersampling */
//		initVd(&vd, traj.FOV, traj.naxes, 0);
//		if(isVD)
//		{
//			addVdStep(&vd, vdType, traj.vdStart, vdParam, 1);
//			addVdStep(&vd, vdfPOLY, traj.vdEnd, 1, traj.FOV[3]/traj.FOV[0]);
//		}
//	}
	
  outputs[0] = mxCreateStructMatrix(1, 1, trajectoryStructFields.size(), trajectoryStructFields.data());
  mxArray* trajectoryStruct = outputs[0];

  mxArray* interpolation = NULL;
  if(interpolationStructFields.size())
    interpolation = mxCreateStructMatrix(1, 1, interpolationStructFields.size(), interpolationStructFields.data());

  if(!generator.generate())
    mexErrMsgTxt("Failed to generate trajectory");

  Trajectory *trajectory = generator.trajectory();
  switch(generator.trajectoryType())
	{
    case CONES:
      {
        Cones* cones = generator.getConesTrajectory();
        copyToStruct(trajectoryStruct, "basisReadoutLengths", cones->numBasisReadoutPoints, trajectory->bases, mxINT32_CLASS, mxREAL);

        memcpy(mxGetData(mxGetField(outputs[0], 0, "basisReadoutLengths")), cones->numBasisReadoutPoints, trajectory->bases*sizeof(int));

        copyToStruct(interpolation, "theta", cones->interpolation.theta, cones->interpolation.readouts, mxSINGLE_CLASS, mxREAL);
        copyToStruct(interpolation, "xyScale", cones->interpolation.scaleXY, cones->interpolation.readouts, mxSINGLE_CLASS, mxREAL);
        copyToStruct(interpolation, "zScale", cones->interpolation.scaleZ, cones->interpolation.readouts, mxSINGLE_CLASS, mxREAL);
        copyToStruct(interpolation, "phi", cones->interpolation.phi, cones->interpolation.readouts, mxSINGLE_CLASS, mxREAL);
        copyToStruct(interpolation, "basisIndex", cones->interpolation.cone, cones->interpolation.readouts, mxINT32_CLASS, mxREAL);
        copyToStruct(interpolation, "readoutIndex", cones->interpolation.readout, cones->interpolation.readouts, mxINT32_CLASS, mxREAL);
        copyToStruct(interpolation, "interleafOnCone", cones->interpolation.interleafOnCone, cones->interpolation.readouts, mxINT32_CLASS, mxREAL);
        copyToStruct(interpolation, "numInterleavesOnCone", cones->interpolation.interleavesOnCone, cones->interpolation.readouts, mxINT32_CLASS, mxREAL);
      }
      mxSetField(trajectoryStruct, 0, "interpolation", interpolation);
			
			break;
    default:
      sprintf(msg, "%d not supported\n", generator.trajectoryType());
      break;
	}

 mwSize numOutputReadouts = generator.storage()==StoreAll ? trajectory->numReadouts : trajectory->bases;
  mwSize trajectoryDimensions[] = {(mwSize)trajectory->numReadoutPoints, (mwSize)trajectory->numDimensions, numOutputReadouts};
  mxSetField(trajectoryStruct, 0, trajectoryFieldName, mxCreateNumericArray(3, trajectoryDimensions, mxSINGLE_CLASS, mxREAL));
        mxArray* array = mxGetField(trajectoryStruct, 0, trajectoryFieldName);
        void* fieldData = mxGetData(array);
        memcpy(fieldData, trajectory->kSpaceCoordinates, mxGetNumberOfElements(array)*sizeof(float));
	
  mwSize densityCompensationDimensions[] = {(mwSize)trajectory->numReadoutPoints, (mwSize)trajectory->numDimensions};
    mxSetField(trajectoryStruct, 0, densityCompensationFieldName, mxCreateNumericArray(2, densityCompensationDimensions, mxSINGLE_CLASS, mxREAL));
     array = mxGetField(trajectoryStruct, 0, densityCompensationFieldName);
     fieldData = mxGetData(array);
      memcpy(mxGetData(array), trajectory->densityCompensation, mxGetNumberOfElements(array)*sizeof(float));
	
    mwSize gradientDimensions[] = {(mwSize)trajectory->numWaveformPoints, (mwSize)trajectory->numDimensions, numOutputReadouts};
    mxSetField(trajectoryStruct, 0, gradientsFieldName, mxCreateNumericArray(3, gradientDimensions, mxSINGLE_CLASS, mxREAL));
        array = mxGetField(trajectoryStruct, 0, gradientsFieldName);
        memcpy(mxGetData(array), trajectory->gradientWaveforms, mxGetNumberOfElements(array)*sizeof(float));

	
  copyToStruct(trajectoryStruct, imageDimensionsFieldName, trajectory->imageDimensions, generator.numDimensions(), mxINT32_CLASS, mxREAL);
		
//	deleteTrajectory(&traj);

	return;
}

