#include "densitycompensation2D.h"

#include "trajectory.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define MAX_LOBES 50
#define LENGTH_LUT	1024

#define DBG_DCFJ 1

float kernelMaxk;
float kernelLUT[LENGTH_LUT];

__inline double P1D(double k, double F)
{
	double val = 1.0;
	double pfk = M_PI*F*k;

	if(k>1e-6)
		val = sin(pfk)*sin(pfk)/(pfk*pfk);

	return val;
}

/* inline float P(float k, float F) */
__inline double P2D(double k, double F)
/*float P2D(float k, float F)*/
{
	double val = 1.0;

	if(k>1e-6)
	{
/* 		val = pow(F*bessj1(F*PI*k)/(2*k), 2)/(PI*PI*F*F*F*F/16);
		val = F*bessj1(F*PI*k)/(2*k)/(PI*F*F/4);*/
		val = F*j1(F*M_PI*k)/(2*k)/(M_PI*F*F/4);
		val *= val;
}

	return val;
}

__inline double P3D(double k, double F)
/*float P3D(float k, float F)*/
{
	double val;
	double s, c;
	double pfk;
	double normVal = 4.5996e+012;

	if(k<1e-6)
		val = 1;
	else
	{
		pfk = M_PI*F*k;
		s = sin(pfk);
		c = cos(pfk);
		val = (s-pfk*c)/(2*M_PI*M_PI*k*k*k) * (s-pfk*c)/(2*M_PI*M_PI*k*k*k)/normVal;
	}

	return val;
}

void setKernelLUT(double F, double maxk, double (*P)(double,double))
{
	int n;
	double k;

	kernelMaxk = maxk;

	for(n=0; n<LENGTH_LUT; n++)
	{
		k = n*maxk/LENGTH_LUT;
		kernelLUT[n] = P(k,F);
	}

	return;
}

__inline double Plut(double k, double Pr)
{
	(void)(Pr);

	double val;
	int idxk;

	idxk = (int)(k/kernelMaxk*LENGTH_LUT);

	val = kernelLUT[idxk];

	return val;
}

size_t multiToSingleIndex(int* indices, int* size, int dimensions)
{
	size_t index = indices[0];
	size_t lowerDimensionSize = size[0];
	int n;
	for(n=1; n<dimensions; n++)
	{
		index += indices[n]*lowerDimensionSize;
		lowerDimensionSize *= size[n];
	}

	return index;
}

void convolutionDensityCompensation(const struct Trajectory *trajectory, int padfront, int* Nall, int numLobes, int numIter, float *W)
{
	int numCt = 40;

	float kCenter[3];
	float kOther[3];
//	float kNormalizationFactor;

	long nptsTotal;
	long countPtsTotal = 0;
	int n;
	int m, p;
	long nn, mm, pp;
	int mmEnd;
	int nnStart, nnEnd;
	int ppStart, ppEnd;
	int ax;
	int intlc, intlo;
	int nc, no;
	double Pcurrent;
	double kDist, kDistSq;
	double F;
	float Pr;
	float Ks[MAX_LOBES];
	double (*P)(double,double);
	int it;
	int doLUT = 1;

	float *WP;	/* convoluion result */
	int *numNeighborPoints;	/* points within kernel radius */

	short***** ctList;
	int *pointsInCompartment;
	int idxCt, idxCto;
	int idxCtoStart;
	int ct[3] = {0,0,0};
	int numCts[3] = {1,1,1};
	long maxCtLength = 0;
	int maxCtOffset;

	nptsTotal = trajectory->readouts*(trajectory->readoutPoints-padfront);
	float fieldOfViewMax = 0;
	float kMax = 0;
	for(n=0; n<trajectory->dimensions; n++)
	{
		fieldOfViewMax = fmax(trajectory->fieldOfView[n], fieldOfViewMax);
		kMax = fmax(5/trajectory->spatialResolution[n], kMax);
	}

	F = fieldOfViewMax*2*kMax;	/* parameter for kernel width */
	float kNormalizationFactor = .5/kMax;

	switch(trajectory->dimensions)
	{
		case 1:
			for(n=0; n<MAX_LOBES; n++)
				Ks[n] = (n+1)/F;
			P = &P1D;
			break;
		case 2:
			Ks[0] = 1.22/F;
			Ks[1] = 2.23/F;
			Ks[2] = 3.24/F;
			P = &P2D;
			break;
		case 3:
			Ks[0] = 1.43/F;
			Ks[1] = 2.46/F;
			Ks[2] = 3.47/F;
			Ks[3] = 4.48/F;
			P = &P3D;
			break;
	}

	Pr = Ks[numLobes-1];	/* kernel radius */

	if(doLUT)
	{
		setKernelLUT(F, Pr, P);
		P = &Plut;
	}

	maxCtOffset = (int)(Pr*numCt+1);

//	WP = matrix2f(trajectory->readouts, N);
	int trajectoryPoints = trajectory->readouts*trajectory->readoutPoints;
	WP = (float*)malloc(trajectoryPoints*sizeof(float));
//	setValMatrix2f(W, trajectory->readouts, N, 1.0);
	for(n=0; n<trajectoryPoints; n++)
		WP[n] = 1;

//	numNeighborPoints = (int**)matrix2(trajectory->readouts, N, 'i');
	numNeighborPoints = (int*)malloc(trajectoryPoints*sizeof(int));
//	setValMatrix2(numNeighborPoints, trajectory->readouts, N, 0, 'i');
	for(n=0; n<trajectoryPoints; n++)
		numNeighborPoints[n] = 1;

	if(DBG_DCFJ)
		printf("Counting points in compartments\n");
	int totalCompartments = 1;
	for(ax=0; ax<trajectory->dimensions; ax++)
	{
		numCts[ax] = numCt;
		totalCompartments *= numCt;
	}

	if(Nall==NULL)
	{
		Nall = (int*)malloc(trajectory->readouts*sizeof(int));
		for(intlc=0; intlc<trajectory->readouts; intlc++)
			Nall[intlc] = trajectory->readoutPoints;
	}

//	pointsInCompartment = (int***)matrix3(numCts[0], numCts[1], numCts[2], 'i');
	pointsInCompartment = (int*)malloc(totalCompartments*sizeof(int));
//	setValMatrix3((void***)pointsInCompartment, numCts[0], numCts[1], numCts[2], 0, 'i');
	for(n=0; n<totalCompartments; n++)
		pointsInCompartment[n] = 0;

	for(intlc=0; intlc<trajectory->readouts; intlc++)
	{
		for(n=padfront; n<Nall[intlc]; n++)
		{
			float kSpaceCoordinates[3];
			trajectoryCoordinates(n, intlc, trajectory, kSpaceCoordinates, NULL);
			for(ax=0; ax<trajectory->dimensions; ax++)
			{
				ct[ax] = round(kSpaceCoordinates[ax]*kNormalizationFactor*numCts[ax]+numCts[ax]/2.0);
				if((ct[ax]<0) || (ct[ax]>numCts[ax]))
				{
					fprintf(stderr, "Assigned compartment out of range\n");
					fprintf(stderr, "\tAxis %d, readout %d, point %d\n", ax, intlc, n);
					fprintf(stderr, "\tct = (%d,%d,%d)\n", ct[0], ct[1], ct[2]);

					fprintf(stderr, "\tk = (%f", kSpaceCoordinates[0]);
					if(trajectory->dimensions>1)
						fprintf(stderr, ",%f", kSpaceCoordinates[1]);
					if(trajectory->dimensions>2)
						fprintf(stderr, "\%f", kSpaceCoordinates[2]);
					fprintf(stderr, ")\n");
					exit(EXIT_FAILURE);
				}
			}
//			pointsInCompartment[ct[0]][ct[1]][ct[2]]++;
			int compartmentIndex = multiToSingleIndex(ct, numCts, trajectory->dimensions);
			pointsInCompartment[compartmentIndex]++;
			maxCtLength = fmax(pointsInCompartment[compartmentIndex], maxCtLength);
			countPtsTotal++;
		}
	}

	printf("Total acquisition points:\t%ld\n", countPtsTotal);
	printf("Maximum points in a compartment:\t%ld\n", maxCtLength);

	ctList = (short*****)malloc(numCts[0]*sizeof(short****));
	for(m=0; m<numCts[0]; m++)
	{
		ctList[m] = (short****)malloc(numCts[1]*sizeof(short***));
		ct[0] = m;
		for(n=0; n<numCts[1]; n++)
		{
			ctList[m][n] = (short***)malloc(numCts[2]*sizeof(short**));
			ct[1] = n;
			for(p=0; p<numCts[2]; p++)
			{
				ct[2] = p;
				int compartmentIndex = multiToSingleIndex(ct, numCts, trajectory->dimensions);
				if(pointsInCompartment[compartmentIndex])
				{
					ctList[m][n][p] = (short**)malloc(pointsInCompartment[compartmentIndex]*sizeof(short*));
					for(idxCt=0; idxCt<pointsInCompartment[compartmentIndex]; idxCt++)
						ctList[m][n][p][idxCt] = (short*)malloc(2*sizeof(short));
				}
			}
		}
	}

	printf("Assigning points to compartments %dx%dx%d\n",numCts[0], numCts[1], numCts[2]);
//	setValMatrix3((void***)pointsInCompartment, numCts[0], numCts[1], numCts[2], 0, 'i');
	for(n=0; n<totalCompartments; n++)
		pointsInCompartment[n] = 0;

	for(intlc=0; intlc<trajectory->readouts; intlc++)
	{
		for(n=padfront; n<Nall[intlc]; n++)
		{
			float kSpaceCoordinates[3];
			trajectoryCoordinates(n, intlc, trajectory, kSpaceCoordinates, NULL);

			/* Get compartment number */
			for(ax=0; ax<trajectory->dimensions; ax++)
				ct[ax] = floor(kSpaceCoordinates[ax]*kNormalizationFactor*numCts[ax]+numCts[ax]/2.0);

			int compartmentIndex = multiToSingleIndex(ct, numCts, trajectory->dimensions);

			ctList[ct[0]][ct[1]][ct[2]][pointsInCompartment[compartmentIndex]][0] = intlc;
			ctList[ct[0]][ct[1]][ct[2]][pointsInCompartment[compartmentIndex]][1] = n;
			pointsInCompartment[compartmentIndex]++;
		}
	}

	/*writeMatrix3((const void***)pointsInCompartment, numCts[0], numCts[1], numCts[2], 'i', "ctlength");*/

	countPtsTotal = 0;
	int compartmentIndicesOther[3];
	for(it=0; it<numIter; it++)
	{
//		setValMatrix2f(WP, trajectory->readouts, N, 0.0);
		for(n=0; n<trajectoryPoints; n++)
			WP[n] = 0;
		for(m=0; m<numCts[0]; m++)
		{
			ct[0] = m;
			if(trajectory->dimensions==2 && DBG_DCFJ)
				printf("ct[%d]\t%.1f%%\n", m, countPtsTotal/(.01*nptsTotal*numIter));
			for(n=0; n<numCts[1]; n++)
			{
				ct[1] = n;
				if(trajectory->dimensions==3 && DBG_DCFJ)
					printf("ct[%d][%d]\t%.1f%%\n", m, n, countPtsTotal/(.01*nptsTotal*numIter));
				for(p=0; p<numCts[2]; p++)
				{
					ct[2] = p;
					int compartmentIndex = multiToSingleIndex(ct, numCts, trajectory->dimensions);
					for(idxCt=0; idxCt<pointsInCompartment[compartmentIndex]; idxCt++)
					{
						countPtsTotal++;

						nc = ctList[m][n][p][idxCt][1];
						intlc = ctList[m][n][p][idxCt][0];
						size_t indexCurrent = intlc*trajectory->readoutPoints + nc;

						float kSpaceCoordinates[3];
						trajectoryCoordinates(nc, intlc, trajectory, kSpaceCoordinates, NULL);

						for(ax=0; ax<trajectory->dimensions; ax++)
							kCenter[ax] = kSpaceCoordinates[ax];

						mmEnd = fmin(m+maxCtOffset+1, numCts[0]);
						for(mm=m; mm<mmEnd; mm++)
						{
							compartmentIndicesOther[0] = mm;
							if(mm==m)
								nnStart = n;
							else
								nnStart = fmax(0,n-maxCtOffset);
							nnEnd = fmin(n+maxCtOffset+1, numCts[1]);
							for(nn=nnStart; nn<nnEnd; nn++)
							{
								compartmentIndicesOther[1] = nn;
								if((nn==n) && (mm==m))
									ppStart = p;
								else
									ppStart = fmax(0, p-maxCtOffset);
								ppEnd = fmin(p+maxCtOffset+1, numCts[2]);
								for(pp=ppStart; pp<ppEnd; pp++)
								{
									compartmentIndicesOther[2] = pp;
									if((mm==m) && (nn==n) && (pp==p))
										idxCtoStart = idxCt;
									else
										idxCtoStart = 0;
									size_t compartmentIndexOther = multiToSingleIndex(compartmentIndicesOther, numCts, trajectory->dimensions);
									for(idxCto=idxCtoStart; idxCto<pointsInCompartment[compartmentIndexOther]; idxCto++)
									{
										no = ctList[mm][nn][pp][idxCto][1];
										intlo = ctList[mm][nn][pp][idxCto][0];
										size_t indexOther = intlo*trajectory->readoutPoints + no;

										float kSpaceCoordinates[3];
										trajectoryCoordinates(no, intlo, trajectory, kSpaceCoordinates, NULL);

										kDistSq = 0;
										for(ax=0; ax<trajectory->dimensions; ax++)
										{
											kOther[ax] = kSpaceCoordinates[ax];
											kDistSq += (kOther[ax]-kCenter[ax])*(kOther[ax]-kCenter[ax]);
										}

										if((no==nc) && (intlo==intlc))
										{
											/* WP[intlc][nc] += PI*PI*F*F*F*F/16*W[intlc][nc];
											//WP[intlc][nc] += 4.5986019e+012*W[intlc][nc]; */
											WP[indexCurrent] += W[indexCurrent];
											/*WP[intlc][nc] += 1.0f;*/
											numNeighborPoints[indexCurrent]++;
										}
										else
										{
											if(kDistSq<(Pr*Pr))
											{
												kDist = sqrt(kDistSq);
												numNeighborPoints[indexCurrent]++;
												numNeighborPoints[indexOther]++;
												/*if(kDist<4e-5)
													//Pcurrent = 4.5986019e+012;
													Pcurrent = 1;
												else*/
													Pcurrent = P(kDist, F);
												/* printf("[%d][%d][%d][%d/%d] - [%d][%d][%d][%d/%d]:\t(%d,%d) - (%d,%d):\t%f\n", m,n,p,idxCt+1,pointsInCompartment[m][n][p], mm,nn,pp,idxCto+1,pointsInCompartment[mm][nn][pp],nc,intlc,no,intlo,Pcurrent);*/
												WP[indexCurrent] += Pcurrent*W[indexOther];
												WP[indexOther] += Pcurrent*W[indexCurrent];
												/*WP[intlc][nc] += 1.0f;
												WP[intlo][no] += 1.0f;*/
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
//		dividePtMatrix2f(W, (const float**)WP, trajectory->readouts, N);
		for(n=0; n<trajectoryPoints; n++)
			W[n] /= WP[n];
	}

	return;
}
