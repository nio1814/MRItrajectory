#include "densitycompensation2D.h"

#include "trajectory.h"

#include <stdlib.h>
#include <math.h>

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

void jdcf(const float ***ktraj, int N, int padfront, int* Nall, int numIntl, int numAxes, int numLobes, float FOV, float kmax, int numCt, int numIter, float **W)
{
	float kCenter[3];
	float kOther[3];
	float kNormFactor;
	float maxw;

	long nptsTotal;
	long countPtsTotal = 0;
	long n, m, p;
	long nn, mm, pp;
	int mmStart, mmEnd;
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

	float **WP;	/* convoluion result */
	int **nbPts;	/* points within kernel radius */

	short***** ctList;
	int ***ctListLength;
	int idxCt, idxCto;
	int idxCtoStart;
	int ct[3] = {0,0,0};
	int numCts[3] = {1,1,1};
	long maxCtLength = 0;
	int maxCtOffset;

	char filename[128];

	nptsTotal = numIntl*(N-padfront);
	F = FOV*2*kmax;	/* parameter for kernel width */
	kNormFactor = 1;	/* assuming data already normalized */

	switch(numAxes)
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

	WP = matrix2f(numIntl, N);
	setValMatrix2f(W, numIntl, N, 1.0);

	nbPts = (int**)matrix2(numIntl, N, 'i');
	setValMatrix2(nbPts, numIntl, N, 0, 'i');

	if(DBG_DCFJ) printf("Counting points in compartments\n");
	for(ax=0; ax<numAxes; ax++)
		numCts[ax] = numCt;

	if(Nall==NULL)
	{
		Nall = (int*)malloc(numIntl*sizeof(int));
		for(intlc=0; intlc<numIntl; intlc++)
			Nall[intlc] = N;
	}

	ctListLength = (int***)matrix3(numCts[0], numCts[1], numCts[2], 'i');
	setValMatrix3((void***)ctListLength, numCts[0], numCts[1], numCts[2], 0, 'i');
	for(intlc=0; intlc<numIntl; intlc++)
	{
		for(n=padfront; n<Nall[intlc]; n++)
		{
			for(ax=0; ax<numAxes; ax++)
			{
				ct[ax] = floor(ktraj[ax][intlc][n]*kNormFactor*numCts[ax]+numCts[ax]/2.0);
				if((ct[ax]<0) || (ct[ax]>numCts[ax]))
				{
					fprintf(stderr, "Assigned compartment out of range\n");
					fprintf(stderr, "\tAxis %d, Cone %d, point %d\n", ax, intlc, n);
					fprintf(stderr, "\tct = (%d,%d,%d)\n", ct[0], ct[1], ct[2]);

					fprintf(stderr, "\tk = (%f", ktraj[0][intlc][n]);
					if(numAxes>1)
						fprintf(stderr, ",%f", ktraj[1][intlc][n]);
					if(numAxes>2)
						fprintf(stderr, "\%f", ktraj[2][intlc][n]);
					fprintf(stderr, ")\n");
					exit(EXIT_FAILURE);
				}
			}
			ctListLength[ct[0]][ct[1]][ct[2]]++;
			maxCtLength = max(ctListLength[ct[0]][ct[1]][ct[2]], maxCtLength);
			countPtsTotal++;
		}
	}

	printf("Total acquisition points:\t%d\n", countPtsTotal);
	printf("Maximum points in a compartment:\t%d\n", maxCtLength);

	ctList = (short*****)malloc(numCts[0]*sizeof(short****));
	for(m=0; m<numCts[0]; m++)
	{
		ctList[m] = (short****)malloc(numCts[1]*sizeof(short***));
		for(n=0; n<numCts[1]; n++)
		{
			ctList[m][n] = (short***)malloc(numCts[2]*sizeof(short**));
			for(p=0; p<numCts[2]; p++)
			{
				if(ctListLength[m][n][p])
				{
					ctList[m][n][p] = (short**)malloc(ctListLength[m][n][p]*sizeof(short*));
					for(idxCt=0; idxCt<ctListLength[m][n][p]; idxCt++)
						ctList[m][n][p][idxCt] = (short*)malloc(2*sizeof(short));
				}
			}
		}
	}

	printf("Assigning points to compartments %dx%dx%d\n",numCts[0],numCts[1],numCts[2]);
	setValMatrix3((void***)ctListLength, numCts[0], numCts[1], numCts[2], 0, 'i');
	for(intlc=0; intlc<numIntl; intlc++)
	{
		for(n=padfront; n<Nall[intlc]; n++)
		{
			/* Get compartment number */
			for(ax=0; ax<numAxes; ax++)
				ct[ax] = floor(ktraj[ax][intlc][n]*kNormFactor*numCts[ax]+numCts[ax]/2.0);

			ctList[ct[0]][ct[1]][ct[2]][ctListLength[ct[0]][ct[1]][ct[2]]][0] = intlc;
			ctList[ct[0]][ct[1]][ct[2]][ctListLength[ct[0]][ct[1]][ct[2]]][1] = n;
			ctListLength[ct[0]][ct[1]][ct[2]]++;
		}
	}

	/*writeMatrix3((const void***)ctListLength, numCts[0], numCts[1], numCts[2], 'i', "ctlength");*/

	countPtsTotal = 0;
	for(it=0; it<numIter; it++)
	{
		setValMatrix2f(WP, numIntl, N, 0.0);
		for(m=0; m<numCts[0]; m++)
		{
			if(numAxes==2 && DBG_DCFJ)
				printf("ct[%d]\t%.1f%%\n", m, countPtsTotal/(.01*nptsTotal*numIter));
			for(n=0; n<numCts[1]; n++)
			{
				if(numAxes==3 && DBG_DCFJ)
					printf("ct[%d][%d]\t%.1f%%\n", m, n, countPtsTotal/(.01*nptsTotal*numIter));
				for(p=0; p<numCts[2]; p++)
				{
					for(idxCt=0; idxCt<ctListLength[m][n][p]; idxCt++)
					{
						countPtsTotal++;

						nc = ctList[m][n][p][idxCt][1];
						intlc = ctList[m][n][p][idxCt][0];

						for(ax=0; ax<numAxes; ax++)
							kCenter[ax] = ktraj[ax][intlc][nc];

						mmEnd = min(m+maxCtOffset+1, numCts[0]);
						for(mm=m; mm<mmEnd; mm++)
						{
							if(mm==m)
								nnStart = n;
							else
								nnStart = max(0,n-maxCtOffset);
							nnEnd = min(n+maxCtOffset+1, numCts[1]);
							for(nn=nnStart; nn<nnEnd; nn++)
							{
								if((nn==n) && (mm==m))
									ppStart = p;
								else
									ppStart = max(0, p-maxCtOffset);
								ppEnd = min(p+maxCtOffset+1, numCts[2]);
								for(pp=ppStart; pp<ppEnd; pp++)
								{
									if((mm==m) && (nn==n) && (pp==p))
										idxCtoStart = idxCt;
									else
										idxCtoStart = 0;
									for(idxCto=idxCtoStart; idxCto<ctListLength[mm][nn][pp]; idxCto++)
									{
										no = ctList[mm][nn][pp][idxCto][1];
										intlo = ctList[mm][nn][pp][idxCto][0];

										kDistSq = 0;
										for(ax=0; ax<numAxes; ax++)
										{
											kOther[ax] = ktraj[ax][intlo][no];
											kDistSq += (kOther[ax]-kCenter[ax])*(kOther[ax]-kCenter[ax]);
										}

										if((no==nc) && (intlo==intlc))
										{
											/* WP[intlc][nc] += PI*PI*F*F*F*F/16*W[intlc][nc];
											//WP[intlc][nc] += 4.5986019e+012*W[intlc][nc]; */
											WP[intlc][nc] += W[intlc][nc];
											/*WP[intlc][nc] += 1.0f;*/
											nbPts[intlc][nc]++;
										}
										else
										{
											if(kDistSq<(Pr*Pr))
											{
												kDist = sqrt(kDistSq);
												nbPts[intlc][nc]++;
												nbPts[intlo][no]++;
												/*if(kDist<4e-5)
													//Pcurrent = 4.5986019e+012;
													Pcurrent = 1;
												else*/
													Pcurrent = P(kDist, F);
												/* printf("[%d][%d][%d][%d/%d] - [%d][%d][%d][%d/%d]:\t(%d,%d) - (%d,%d):\t%f\n", m,n,p,idxCt+1,ctListLength[m][n][p], mm,nn,pp,idxCto+1,ctListLength[mm][nn][pp],nc,intlc,no,intlo,Pcurrent);*/
												WP[intlc][nc] += Pcurrent*W[intlo][no];
												WP[intlo][no] += Pcurrent*W[intlc][nc];
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
		dividePtMatrix2f(W, (const float**)WP, numIntl, N);
		/*sprintf(filename, "data/W%d", it+1);
		writeMatrix2f((const float**)W, numIntl, N, filename);
		sprintf(filename, "data/WP%d", it+1);
		writeMatrix2f((const float**)WP, numIntl, N, filename);
		/*maxw = maxVal2f((const float**)W, numIntl, N);
		scaleMatrix2f(W, numIntl, N, 1/maxw);
		maxw = maxVal2f((const float**)W, numIntl, N);*/
	}

	/*writeMatrix3((const void***)ctListLength, numCt, numCt, numCts[2], 'i', "ctlistlength");
	writeMatrix2((const void**)nbPts, numIntl, N, 'i', "nbpts");*/

	return;
}
