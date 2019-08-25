#include "spinwarp.h"

//#include <string.h>
//#ifndef HW_IO
//#include <stdio.h>
//#endif

#include "variabledensity.h"
#include "trajectory.h"
#include "mrgradient.h"
#include "mathops.h"
#include "phaseencode.h"
#include "arrayops.h"
//#include "mtg_functions.h"
//#include "trajectory.h"
//#include "vdrad.h"
//#include "cmatrix.h"
//#include "dataformat.h"

//#ifndef max
//#define max(a,b) ((a) < (b) ? (b) : (a))
//#endif

//#ifndef min
//#define min(a,b) ((a) < (b) ? (a) : (b))
//#endif

//#include <math.h>
#include <stdlib.h>
#include <string.h>


struct CartesianSchedule* newCartesianSchedule(const int length, const int numPhaseEncodeAxes)
{
  if(length < 1)
  {
    fprintf(stderr, "newCartesianSchedule: schedule length %d < 1", length);
    return NULL;
  }

  struct CartesianSchedule* schedule = (struct CartesianSchedule*)malloc(sizeof(struct CartesianSchedule));
  schedule->phaseEncodeIndex1 = (int*)malloc(length*sizeof(int));
  if(numPhaseEncodeAxes>1)
    schedule->phaseEncodeIndex2 = (int*)malloc(length*sizeof(int));
  schedule->length = length;

  return schedule;
}

//    return status;
//}

//int spinwarpLoadInfo(const char* filename, struct SpinwarpInfo *info)
//{
//    int status = 1;
//    FILE* file;
//    int ver = 0;
//    int doSwapE = doSwapEndian('l');
//    int N;

//    file = fopen(filename, "rb");
//    if(file==NULL)
//	   fprintf(stderr,"spinwarpSavInfo: Error opening file %s\n", filename);
//    else
//    {
//	   fread(&ver, sizeof(int), 1, file);
//	   if(doSwapE)
//		   swap_int_buffer(&ver, 1);

//	   fread(&info->ndim, sizeof(int), 1, file);
//	   fread(info->npe, sizeof(int), 2, file);
//	   if(doSwapE)
//	   {
//		   swap_int_buffer(&info->ndim, 1);
//		   swap_int_buffer(info->npe, 2);
//	   }

//	   if(info->ndim>2)
//		N = info->npe[0]*info->npe[1]*(info->ndim-1);
//	else
//		N = info->npe[0]*(info->ndim-1);

//	   info->pes = (float*)malloc(N*sizeof(float));
//	   fread(info->pes, sizeof(float), N, file);
//	   fread(info->fov, sizeof(float), 6, file);
//	   fread(info->res, sizeof(float), 3, file);
//	   fread(&info->ndep, sizeof(int), 1, file);
//	   fread(&info->nread, sizeof(int), 1, file);
//	   fread(&info->nrep, sizeof(int), 1, file);
//	   fread(info->axid, sizeof(int), 3, file);
//	   fread(&info->dt, sizeof(float), 1, file);

//		if(doSwapE)
//		{
//			swap_float_buffer(info->pes, N);
//			swap_float_buffer(info->fov, 6);
//			swap_float_buffer(info->res, 3);
//		   swap_int_buffer(&info->ndep, 1);
//		   swap_int_buffer(&info->nread, 1);
//		   swap_int_buffer(&info->nrep, 1);
//		   swap_int_buffer(info->axid, 3);
//			swap_float_buffer(&info->dt, 1);
//		}

//		   N = (info->ndep+info->nread+info->nrep)*info->ndim;
//		   info->g = (float*)malloc(N*sizeof(float));
//		  fread(info->g, sizeof(float), N, file);

//		if(doSwapE)
//			swap_float_buffer(info->g, N);

//	   fclose(file);
//	   status = 0;
//    }

//    return status;
//}

//int spinwarpSaveSched(const char* filename, struct SpinwarpSched *sched)
//{
//	int status = 1;
//	FILE *file;
//	int doSwapE;
//    int ver = 0;
//	int s;

//	file = fopen(filename, "wb");
//	if(file!=NULL)
//	{
//		doSwapE = doSwapEndian('l');

//		if(doSwapE)
//		{
//		  swap_int_buffer(&ver, 1);
//			swap_int_buffer(&sched->ndim, 1);
//			swap_int_buffer(&sched->length, 1);
//		}

//	   fwrite(&ver, 1, sizeof(int), file);

//		fwrite(&sched->ndim, 1, sizeof(int), file);
//		fwrite(&sched->length, 1, sizeof(int), file);

//		if(doSwapE)
//		{
//			swap_int_buffer(&sched->ndim, 1);
//			swap_int_buffer(&sched->length, 1);

//			swap_int_buffer(sched->instruction, sched->length);
//			swap_float_buffer(sched->s1, sched->length);
//			if(sched->ndim>2)
//				swap_float_buffer(sched->s2, sched->length);
//		}

//		for(s=0; s<sched->length; s++)
//		{
//			fwrite(&sched->instruction[s], sizeof(int), 1, file);
//			fwrite(&sched->s1[s], sizeof(float), 1, file);
//			if(sched->ndim>2)
//				fwrite(&sched->s2[s], sizeof(float), 1, file);
//		}

//		if(doSwapE)
//		{
//			swap_int_buffer(sched->instruction, sched->length);
//			swap_float_buffer(sched->s1, sched->length);
//			if(sched->ndim>2)
//				swap_float_buffer(sched->s2, sched->length);
//		}

//		fclose(file);

//		status = 0;
//	}
//	else
//		fprintf(stderr, "radial_save_sched() Error loading %s\n", filename);

//	return status;
//}

//int spinwarpLoadSched(const char* filename, struct SpinwarpSched *sched)
//{
//	int status = 1;
//	int doSwapE = doSwapEndian('l');
//	FILE *file;
//	int ver = 0;
//	int s;

//	file = fopen(filename, "rb");
//	if(file!=NULL)
//	{
//		fread(&ver, sizeof(int), 1, file);

//		fread(&sched->ndim, sizeof(int), 1, file);
//		fread(&sched->length, sizeof(int), 1, file);

//		if(doSwapE)
//		{
//		  swap_int_buffer(&ver, 1);
//			swap_int_buffer(&sched->ndim, 1);
//			swap_int_buffer(&sched->length, 1);
//		}

//		initNDFTsched(sched, sched->length, sched->ndim);

//		for(s=0; s<sched->length; s++)
//		{
//			fread(&sched->instruction[s], sizeof(int), 1, file);
//			fread(&sched->s1[s], sizeof(float), 1, file);
//			if(sched->ndim>2)
//				fread(&sched->s2[s], sizeof(float), 1, file);
//		}

//		if(doSwapE)
//		{
//			swap_int_buffer(sched->instruction, sched->length);
//			swap_float_buffer(sched->s1, sched->length);
//			if(sched->ndim>2)
//				swap_float_buffer(sched->s2, sched->length);
//		}

//		fclose(file);
//		status = 0;
//	}
//	else
//		fprintf(stderr, "radial_save_sched() Error loading %s\n", filename);

//	return status;
//}

void cartesianReadoutGradients(const int readout, const struct CartesianSchedule *schedule, const struct Cartesian* cartesian, float *gx, float *gy, float *gz)
{
  float *gradientWaveform[3] = {gx, gy, gz};

  memcpy(&gradientWaveform[cartesian->readoutAxis][cartesian->numPreReadoutWaveformPoints], cartesian->readoutGradient, cartesian->numReadoutWaveformPoints*sizeof(float));
  for(int p=0; p<cartesian->numPhaseEncodeAxes; p++)
  {
    const struct PhaseEncoder* phaseEncoder = cartesian->phaseEncoders[p];
    int axis = cartesian->phaseEncodeAxes[p];
    int offset = cartesian->numPrePhaseEncodeWaveformPoints[p];
    memcpy(&gradientWaveform[axis][offset], phaseEncoder->gradient, phaseEncoder->numPoints*sizeof(float));

    const int* phaseEncodeIndices = p ? schedule->phaseEncodeIndex2 : schedule->phaseEncodeIndex1;
    const int phaseEncodeIndex = phaseEncodeIndices[readout];

    const float scale = phaseEncoder->scales[phaseEncodeIndex];
    scalefloats(gradientWaveform[axis], phaseEncoder->numPoints, scale);
  }
}

/**
 \param[in]	ndep	number of points before readout begins including phase encode and readout dephase and ramp up
 \param[in]	ndep	number of points after readout end including phase encode and readout ramp down
*/
struct Cartesian generateCartesian(float* fieldOfView, float *spatialResolution, const char* axisOrder, int numDimensions, int doRephase, float samplingInterval, float gradientLimit, float slewLimit)
{
  const int readoutAxis = axisNameToIndex(axisOrder[0]);
  fieldOfView[readoutAxis] = MAX(calculateMinFieldOfView(gradientLimit, samplingInterval), fieldOfView[readoutAxis]);

  int imageDimensions[3];
  for(int d=0; d<numDimensions; d++)
    adjustSpatialResolution(fieldOfView[d], &imageDimensions[d], &spatialResolution[d]);

  struct Cartesian cartesian;
  cartesian.samplingInterval = samplingInterval;
  cartesian.readoutSpatialResolution = spatialResolution[readoutAxis];

	/* make readout gradient */
  cartesian.numReadoutPoints = imageDimensions[readoutAxis];
  int numReadoutRampPoints;	/* number of ramp points for readout gradient */
  int numReadoutDephasePoints;	/* number of points before readout begins for readout gradient */
  grd_readout(1, doRephase, cartesian.numReadoutPoints, spatialResolution[readoutAxis], gradientLimit, slewLimit, samplingInterval, &cartesian.readoutGradient, &numReadoutRampPoints, &numReadoutDephasePoints, &cartesian.numReadoutWaveformPoints);

  cartesian.numPhaseEncodeAxes = numDimensions - 1;

  int numPhaseEncodePoints = 0;
  cartesian.numWaveformPoints = cartesian.numReadoutWaveformPoints;
  for(int p=0; p<cartesian.numPhaseEncodeAxes; p++)
  {
    const int axis = axisNameToIndex(axisOrder[p+1]);
    cartesian.phaseEncodeAxes[p] = axis;
    cartesian.phaseEncoders[p] = newPhaseEncoder(fieldOfView[axis], spatialResolution[axis], 1, 1, doRephase, cartesian.numReadoutPoints, samplingInterval, gradientLimit, slewLimit);
    numPhaseEncodePoints = MAX(numPhaseEncodePoints, cartesian.phaseEncoders[p]->numPointsPhaseEncode);

    cartesian.numWaveformPoints = MAX(cartesian.phaseEncoders[p]->numPoints, cartesian.numWaveformPoints);
    cartesian.phaseEncodeFieldOfView[p] = fieldOfView[axis];
    cartesian.phaseEncodeSpatialResolution[p] = spatialResolution[axis];
  }

  cartesian.numPreReadoutPoints = MAX(numReadoutDephasePoints, numPhaseEncodePoints);
  cartesian.numPreReadoutWaveformPoints = cartesian.numPreReadoutPoints - numReadoutDephasePoints;
  for(int p=0; p<cartesian.numPhaseEncodeAxes; p++)
  {
    cartesian.numPrePhaseEncodeWaveformPoints[p] = cartesian.numPreReadoutPoints- cartesian.phaseEncoders[p]->numPointsPhaseEncode;
  }
//  *numRephasePoints = doRephase ? MAX(numReadoutDephasePoints, numPhaseEncodePoints) : 0;

//  cartesian.numWaveformPoints = *numDephasePoints + cartesian.numReadoutWaveformPoints + *numRephasePoints;

//	*g = (float*)calloc(ndim*npts, sizeof(float));
//  TODO make phase encodes and ramp time the same
//  if(numReadoutRampPoints > numPhaseEncodePoints)
//	{
//		/* readout dephase + ramp up longer than phase encode
//			recreate longer phase encode with reduced amplitude */
//    const float scale = MIN(numPhaseEncodePoints/(float)(numReadoutRampPoints), 1.0f);

//		gmax2 = scale*gp[npr];
//		smax2 = scale*sp[1];

//		free(gp);
//		free(sp);
//    grd_trapezoid(gmax2, smax2, kmaxpe, dt, NULL, &gp, &sp, NULL, &np, &npr);

//		memcpy(*g, gr, npts*sizeof(float));
//		for(n=1; n<ndim; n++)
//			memcpy(&(*g)[n*npts+ndr-np], gp, np*sizeof(float));
//	}
//	else
//	{
//		memcpy(&((*g)[np-ndr]), gr, (ndr+*nread+ndr)*sizeof(float));
//		for(n=1; n<ndim; n++)
//			memcpy(&(*g)[n*npts], gp, np*sizeof(float));
//	}

	/* create second phase encode gradient to move back to k=0 */
//	for(n=1; n<ndim; n++)
//	{
//		ptr = &(*g)[n*npts+*ndep+*nread];
//		memcpy(ptr, gp, np*sizeof(float));
//		scalefloats(ptr, np, -1.0f);	/* flip polarity */
//	}

//	switch(peOrder)
//	{
//		case swRASTER:
//			if(vd->nsteps==1)
//			{
//				npes[0] = imgN[axid[1]];
//				p1 = (float*)malloc(npes[0]*sizeof(float));
//				calcPeScales(npes[0], p1, 1, 1);
//				if(kmaxs[axid[1]]<kmaxpe)
//					scalefloats(p1, npes[0], kmaxs[axid[1]]/kmaxpe);

//				if(ndim>2)
//				{
//					npes[1] = imgN[axid[2]];
//					p2 = (float*)malloc(npes[1]*sizeof(float));
//					calcPeScales(npes[1], p2, 1, 1);
//					if(kmaxs[axid[2]]<kmaxpe)
//						scalefloats(p2, npes[1], kmaxs[axid[2]]/kmaxpe);
//				}
//				else
//					npes[1] = 1;

//				*npe = npes[0]*npes[1];
//				*pes = (float*)malloc(*npe*(ndim-1)*sizeof(float));

//				for(m1=0; m1<npes[0]; m1++)
//					for(m2=0; m2<npes[1]; m2++)
//					{
//						n = m2*npes[0] + m1;
//						if(ndim>2)
//						{
//							(*pes)[2*n] = p1[m1];
//							(*pes)[2*n+1] = p2[m2];
//						}
//						else
//							(*pes)[n] = p1[m1];
//					}

//				if(doCornerCut)
//				{

//				}
//			}
//			else
//			{
//			}
//			break;
//		case swVDRAD:
//			/*npes[0] = imgN[axid[1]]/res[axid[1]];
//			npes[1] = imgN[axid[2]]/res[axid[2]];*/
//			npes[0] = imgN[axid[1]];
//			npes[1] = imgN[axid[2]];

//			i1 = (int*)malloc(npes[0]*npes[1]*sizeof(int));
//			i2 = (int*)malloc(npes[0]*npes[1]*sizeof(int));

//			scale = sqrt(1/vd->scale[vd->nsteps-1]);

//			/* get calibratoin region size */
//			if(vd->nsteps==1)
//				n = 0;
//			else
//				/*if(kmaxs[axid[1]]>kmaxs[axid[2]])
//					n = vd->kr[1]/kmaxpe*imgN[axid[1]];
//				else
//					n = vd->kr[1]/kmaxpe*imgN[axid[2]];*/

//				n = 0;
//				for(m1=1; m1<=2; m1++)
//					n = fmaxf(n, vd->kr[1]/kmaxs[axid[m1]]*imgN[axid[m1]]);

//			if(nPerSeg<=0)
//				nPerSeg = fmaxf(imgN[axid[1]], axid[2])/4;

//			genRadialSamplingOrdering(i1, i2, 10*vd->fov[axid[1]], 10*vd->fov[axid[2]], npes[0], npes[1], scale, scale, sqrt(2.0)*n, nPerSeg, doCornerCut, npes[0]*npes[1], 0, 40, 1, npe);

//			*pes = (float*)malloc(*npe*2*sizeof(float));
//			for(n=0; n<*npe; n++)
//			{
//				(*pes)[2*n] = (2.0f*i1[n]-imgN[axid[1]])/imgN[axid[1]];
//				(*pes)[2*n+1] = (2.0f*i2[n]-imgN[axid[2]])/imgN[axid[2]];
//			}

//			npes[0] = *npe;
//			npes[1] = 1;

//			break;
//	}

//	if(p1)
//		free(p1);
//	if(p2)
//		free(p2);
//	free(gr);
//	free(sp);
//	free(gp);
//	if(i1)
//		free(i1);
//	if(i2)
//		free(i2);

  return cartesian;
}

struct Trajectory* cartesianToTrajectory(const struct Cartesian* cartesian, const enum WaveformStorageType storage)
{
//	if(vd->nsteps==1)
//		peOrder = swRASTER;

  int numReadouts = 1;
  for(int p=0; p<cartesian->numPhaseEncodeAxes; p++)
    numReadouts *= cartesian->phaseEncoders[p]->numPhaseEncodes;

  struct Trajectory* trajectory = newTrajectory();

  allocateTrajectory(trajectory, cartesian->numReadoutPoints, cartesian->numWaveformPoints, 1 + cartesian->numPhaseEncodeAxes, 1, numReadouts, storage);
  trajectory->fieldOfView[cartesian->readoutAxis] = cartesian->readoutFieldOfView;
  trajectory->spatialResolution[cartesian->readoutAxis] = cartesian->readoutSpatialResolution;
  for(int p=0; p<cartesian->numPhaseEncodeAxes; p++)
  {
    trajectory->fieldOfView[cartesian->phaseEncodeAxes[p]] = cartesian->phaseEncodeFieldOfView[p];
    trajectory->spatialResolution[cartesian->phaseEncodeAxes[p]] = cartesian->phaseEncodeSpatialResolution[p];
  }

  trajectory->samplingInterval = cartesian->samplingInterval;

  struct CartesianSchedule* schedule = newCartesianSchedule(trajectory->numReadouts, cartesian->numPhaseEncodeAxes);
  for(int n=0; n<schedule->length; n++)
	{
    schedule->phaseEncodeIndex1[n] = n % cartesian->phaseEncoders[0]->numPhaseEncodes;
    if(schedule->phaseEncodeIndex2)
      schedule->phaseEncodeIndex2[n] = n / cartesian->phaseEncoders[0]->numPhaseEncodes;
	}

  if(storage==STORE_ALL)
	{
    setValues(trajectory->densityCompensation, trajectory->numReadoutPoints*trajectory->numReadouts, 1.0f);

    for(int r=0; r<trajectory->numReadouts; r++)
		{
      float* gradientWaveformX = trajectoryGradientWaveform(trajectory, r, 0);
      float* gradientWaveformY = trajectoryGradientWaveform(trajectory, r, 1);
      float* gradientWaveformZ = cartesian->numPhaseEncodeAxes > 1 ? trajectoryGradientWaveform(trajectory, r, 2) : NULL;
      cartesianReadoutGradients(r, schedule, cartesian, gradientWaveformX, gradientWaveformY, gradientWaveformZ);
		}
	}
	else
	{
    float* gradientWaveform = trajectoryGradientWaveform(trajectory, 0, cartesian->readoutAxis);
    memcpy(gradientWaveform, cartesian->readoutGradient, cartesian->numReadoutWaveformPoints*sizeof(float));
    for(int p=0; p<cartesian->numPhaseEncodeAxes; p++)
    {
      const int axis = cartesian->phaseEncodeAxes[p];
      gradientWaveform = trajectoryGradientWaveform(trajectory, 0, axis);
      memcpy(gradientWaveform, cartesian->phaseEncoders[p]->gradient, cartesian->phaseEncoders[p]->numPoints * sizeof(float));
    }
	}

  return trajectory;
}

