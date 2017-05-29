#include "spiral.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "trajectory_s.h"
#include "cmatrix.h"
#include "mathfunctions.h"
#include "spiralgen.h"
#include "bh_array.h"

#ifndef MAX_PG_IAMP
#define MAX_PG_IAMP 32766
#endif

#define DBG_SPR 1

#ifdef MGD_TGT
#include "mgdmath.h"
#endif

void makeSpiralInterleavesTraj(struct Trajectory* spiral, int doK, int doG, int doI)
{
	int intl, i;
	float phi, cosphi, sinphi;

	for(intl=0; intl<spiral->ninter; intl++)
	{
	     /* for XY rotation of intl:0~Nintl-1 */
	     phi = 2*PI*(float)intl/(float)spiral->ninter;
	     cosphi = cos( phi );
	     sinphi = sin( phi );

         for(i=0; i<spiral->npts; i++)
	     {
			if(doK && (i<spiral->acqLength))
			{
				spiral->ks[0][intl][i] = cosphi*spiral->ks[0][0][i] - sinphi*spiral->ks[1][0][i];
				spiral->ks[1][intl][i] = sinphi*spiral->ks[0][0][i] + cosphi*spiral->ks[1][0][i];
			}
			
			if(doG)
			{
				spiral->waveformsG[0][intl][i] = cosphi*spiral->waveformsG[0][0][i] - sinphi*spiral->waveformsG[1][0][i];
				spiral->waveformsG[1][intl][i] = sinphi*spiral->waveformsG[0][0][i] + cosphi*spiral->waveformsG[1][0][i];
			}
			
			if(doI)
			{
				spiral->waveformsI[0][intl][i] = (short)((cosphi*(spiral->waveformsI[0][0][i]) - sinphi*(spiral->waveformsI[1][0][i]))) & 0xfffe;
				spiral->waveformsI[1][intl][i] = (short)((sinphi*(spiral->waveformsI[0][0][i]) + cosphi*(spiral->waveformsI[1][0][i]))) & 0xfffe;
			}
		 }

		 if(doI)
		 {
			/*Set last instruction to odd value*/
			i = spiral->npts-1;
			
            spiral->waveformsI[0][intl][i] |=  0x0001;
            spiral->waveformsI[1][intl][i] |=  0x0001;
		}
	}

	return;
}

void makeInterleaves(float** wavesInit, int ninter, int npts, float*** wavesAll)
{
	int intl, i;
	float phi, cosphi, sinphi;

	for(intl=0; intl<ninter; intl++)
	{
	     /* for XY rotation of intl:0~Nintl-1 */
		 phi = 2*PI*(float)intl/(float)ninter;
	     cosphi = cos(phi);
	     sinphi = sin(phi);

		 for(i=0; i<npts; i++)
	     {
			/* X-grad */
			wavesAll[0][intl][i] = (float)((cosphi*(wavesInit[0][i]) - sinphi*(wavesInit[1][i])));

			/* Y-grad */
			wavesAll[1][intl][i] = (float)((sinphi*(wavesInit[0][i]) + cosphi*(wavesInit[1][i])));
	     }
	}

	return;
}

void makeInterleavesK(struct Trajectory* spiral)
{
	int intl, i;
	float phi, cosphi, sinphi;

	for(intl=0; intl<spiral->ninter; intl++)
	{
	     /* for XY rotation of intl:0~Nintl-1 */
	     phi = 2*PI*(float)intl/(float)spiral->ninter;
	     cosphi = cos( phi );
	     sinphi = sin( phi );

		 for(i=1; i<spiral->acqLength; i++)
	     {
			/* X-grad */
			spiral->ks[0][intl][i] = cosphi*spiral->ks[0][0][i] - sinphi*spiral->ks[1][0][i];
			
			/* Y-grad */
			spiral->ks[1][intl][i] = sinphi*spiral->ks[0][0][i] + cosphi*spiral->ks[1][0][i];
		 }
	}

	return;
}

int spiralLoadTrajRth(char* ksfile, char* wavfile, struct Trajectory *traj, char emodek, char emodew)
{
	int status = 1;
	
/*	Trajectory must be initialized (null memory) or previously allocated or else error will result*/
	deleteTrajectory(traj);
	
	printf("Loading wav file %s\n", wavfile);
	status = loadRthHeader(wavfile, traj, emodew, NULL, 1);
	
	if(!status)
		return status;
		
	traj->ks = (float***)matrix3f(traj->naxes, traj->ninter, traj->acqLength);
	traj->dcf = (float**)matrix2f(traj->ninter, traj->acqLength);
	
	printf("Loading ks file %s\n", ksfile);
	status = loadks(ksfile, traj->ks, traj->dcf, traj->acqLength, traj->ninter, traj->naxes, emodek);
	
	traj->trajType = tjSPIRAL;
	
	return status;
}

/*void initSpiral(struct Trajectory* spiral)
{
	if(spiral->kX==NULL)
		spiral->kX = (float*)malloc(spiral->acqLength*sizeof(float));

	wave2k(spiral->waveformsI[0][0], spiral->waveformsI[1][0], spiral->kX, spiral->kY, spiral->acqLength);

	rewindTraj(&spiral->kX, &spiral->kY, &spiral->GX, &spiral->GY, spiral->max_G, spiral->max_slew, spiral->acqLength, &spiral->npts);
	makeInterleavesG(spiral,0);
	makeInterleavesk(spiral);
	vecdcf(spiral->GX, spiral->GY, spiral->kX, spiral->kY, spiral->acqLength, spiral->kmax, spiral->dcf);

	return;
}*/

void calcSpiralDcf(float *gx, float *gy, float *kx, float *ky, int rolen, float *denscomp)
{
	float *kr, *gr;
	float* Gtwist;
	float *denscomp_it, *denscomp_is;
	float *kPhasexy, *kMagxy;
	int n;

	/* Calculate k-space and gradient magnitudes */
	kr = (float*)malloc(rolen*sizeof(float));
	gr = (float*)malloc(rolen*sizeof(float));
	for(n=0; n<rolen; n++)
	{
        kr[n] = sqrt(kx[n]*kx[n] + ky[n]*ky[n]);
        gr[n] = sqrt(gx[n]*gx[n] + gy[n]*gy[n]);
	}

	Gtwist = (float*)malloc(rolen*sizeof(float));
	kPhasexy = (float*)malloc(rolen*sizeof(float));
	kMagxy = (float*)malloc(rolen*sizeof(float));

	calcphase(kx, ky, kPhasexy, rolen, 0, 1.0);
	calcmag(kx, ky, kMagxy, rolen);
	unwrapPhase(kPhasexy, Gtwist, rolen, rolen, PI);
	for(n=0; n<rolen-1; n++)
        Gtwist[n] = max((Gtwist[n+1]-Gtwist[n])/(kr[n+1]-kr[n])*kr[n], 0.0f);
	
	Gtwist[rolen-1] = Gtwist[rolen-2];

	/*% Density compensation due to inter-trajectory spacing (ignoring NINT)
	denscomp_it = abs(kcxy)./sqrt(1+Gtwist.^2);*/
	denscomp_it = (float*)malloc(rolen*sizeof(float));
	for(n=0; n<rolen; n++)
    {
        if(Gtwist[n]>=0)
            denscomp_it[n] = kMagxy[n]/sqrt(1+Gtwist[n]*Gtwist[n]);
        else
            denscomp_it[n] = 0;
    }

	/*% Density compensation due to inter-sample spacing
	denscomp_is = (gr+[gr(2:end); gr(end)])/2;*/
	denscomp_is = (float*)malloc(rolen*sizeof(float));
	addfloats(gr, &(gr[1]), denscomp_is, rolen-1);
	denscomp_is[rolen-1] = 2*gr[rolen-1];
	scalefloats(denscomp_is, rolen, 0.5);

	multiplyfloats(denscomp_is, denscomp_it, denscomp, rolen);

	/*Deallocate mem*/
	free(kr);
	free(gr);
	free(Gtwist);
	free(denscomp_it);
	free(denscomp_is);
	free(kPhasexy);
	free(kMagxy);

	return;
}

float calcc1Spiral(float fov, float nintl)
{
    return 4.0f*PI*PI*fov*fov/(nintl*nintl);
}

float calcGtwistSpiral(float kr, float fr, float nintl)
{
    float gtwist;
    float c1 = calcc1Spiral(fr, nintl);
    /*float c2 = 1.0f;*/	/* only twist when necessary */
    float c2 = 0.0f;	/* always twist */

    gtwist = c1*kr*kr - c2;
    gtwist = sqrt((gtwist>0)?gtwist:0);

    return gtwist;
}

float calcFovFermatFloret(float fov, float kr, float a0)
{
	return 2*a0*fov*fov*kr;
}

int genspiral(float res, struct Vd *vd, int nintl, int maxlen, float dt, enum SpiralType sptype, float floretAngle, float gmax, float smax, float *kstart, int *ntwist, float **kout, float **gout, int *Nk)
{
    int osr = rnd(dt/SPIRALDT);
    float gamrast = GMR*SPIRALDT; /* gamrast*g = dk */
    float dgc = smax*SPIRALDT;	/* the most the gradients can change in 1 raster period */

    float phi = 0;
    float dphi;
    float phiNext;

    float kmax = 5/res;
    float kr;
    float dkr;
    float dkr1, dkr2;
    float dpdk;
    /*float knorm;*/
    int n;
    float *k;
    float g[3];
    int *gsign;
    float gs;
    float gxy;

    float knext[3];
    float dk[3];
    float u[3];
    float ugxy;
    float uxym;
    float gmrange[2];
    float gm;
    float gscale;
    float gmaxread;
    int ndim=2;
    int d;
    float gtwist;
    float term;
    float fov;
    int bt; /* 1 when backtracking */

    int isgood = 1;
    int m;

    maxlen *= osr;
    k = (float*)calloc(ndim*maxlen, sizeof(float));
    gsign = (int*)malloc(maxlen*sizeof(int));

    for(n=0; n<maxlen; n++)
        gsign[n] = 1;

    k[0] = gamrast*dgc;
    k[2] = 2*gamrast*dgc;

    kr = k[2];

    dkr1 = k[0];
    dkr2 = k[2] - k[0];

    n = 1;
    m = 1;
    *ntwist = -1;

    while((kr <= kmax) && isgood)
    {
        /*for(d=0; d<ndim; d++)
            km[d] = 1.5*k[n*ndim+d] - 0.5*k[(n-1)*ndim+d];*/
		  dkr1 = norm2(&k[n*ndim], ndim) - norm2(&k[(n-1)*ndim], ndim);
		  if(n>=2)
			dkr2 = norm2(&k[(n-1)*ndim], ndim) - norm2(&k[(n-2)*ndim], ndim);
		  else
			  dkr2 = 0;
	/*`	dkr = 2*dkr1 - dkr2;*/
		dkr = dkr1;

	   /*kmr = norm2(km, ndim);
	   knorm = kr/kmax;*/
        getFov(vd, kr, &fov);
	   gmaxread = min(maxFilterGrad(fov, dt), gmax);

        if(sptype==spFERMAT)
        {
            fov = calcFovFermatFloret(fov, kr, floretAngle);
            gtwist = 2*PI*kr*fov/nintl;
        }
        else
            gtwist = calcGtwistSpiral(kr, fov, nintl);

        phi = atan2(k[n*ndim+1], k[n*ndim]);
        dpdk = gtwist/kr;
        /*dkr = kmr-kr;*/
        dphi = dpdk*dkr;
        phiNext = phi+dphi;

        knext[0] = cos(phiNext);
        knext[1] = sin(phiNext);
        /*scalefloats(knext, ndim, kmr);*/
        scalefloats(knext, ndim, kr+dkr);
        subtractArray(knext, &k[n*ndim], dk, ndim);
        copyfloats(dk, u, ndim);
        scalefloats(u, ndim, 1.0f/norm2(dk,ndim));

        for(d=0; d<ndim; d++)
            g[d] = (k[n*ndim+d]-k[(n-1)*ndim+d])/gamrast;

        uxym = sqrt(u[0]*u[0]+u[1]*u[1]);
        ugxy = dot(u, g, 2);
        term = uxym*uxym*(dgc*dgc - (g[0]*g[0] + g[1]*g[1])) + ugxy*ugxy;

        if(term>=0)
        {
            for(d=0; d<2; d++)
            {
                if(d==0)
                    gs = -1.0f;
                else
                    gs = 1.0f;

                gmrange[d] = (ugxy + gs*sqrt(term))/(uxym*uxym);
            }

            if(gmrange[0]>gmrange[1])
                bt = 1;
            else
                bt = 0;
        }
        else
            bt = 1;

        if(!bt)
        {
            if(gsign[n]==1)
                gm = gmrange[1];
            else
                gm = gmrange[0];

            gxy = gm*uxym;

            if(gxy>gmax)
                gscale = gmax/gxy;
            else
                gscale = 1.0f;

            gm *= gscale;

            for(d=0; d<ndim; d++)
            {
                g[d] = gm*u[d];
                k[(n+1)*ndim+d] = k[n*ndim+d] + g[d]*gamrast;
            }

            if(gtwist && *ntwist==-1)
                *ntwist = n;

            n++;
        }
        else
        {
            while((n>3) && (gsign[n-1]==-1))
                n--;
            gsign[n-1] = -1;
            n -= 2;
        }
        kr = norm2(&k[n*ndim], ndim);

        isgood = n<(maxlen-1);
        m++;
    }

    if(isgood)
    {
        *Nk = (n-1)/osr;

        *kout = (float*)malloc(*Nk*ndim*sizeof(float));
        *gout = (float*)malloc(*Nk*ndim*sizeof(float));

        for(d=0; d<ndim; d++)
        {
            (*kout)[*Nk*d] = k[d];
            (*gout)[*Nk*d] = k[d]/gamrast;
        }

        for(n=1; n<*Nk; n++)
        {
            for(d=0; d<ndim; d++)
            {
                (*kout)[*Nk*d+n] = k[osr*n*ndim+d];
                (*gout)[*Nk*d+n] = (k[osr*n*ndim+d]-k[osr*(n-1)*ndim+d])/gamrast/osr;
            }
        }
    }

    free(k);
    free(gsign);

    return isgood;
}

void makeSpiralT(struct Vd *vd, float res, float Tdes, float dt, enum SpiralType sptype, float floretAngle, float fovFilt, float gmax, float smax, int *Nk, int *Ng, int *nintl, float *gpre, int nptsGpre, int ndimGpre, int doOutputAll, float **gAll, float **kAll, float **dcfAll)
{
    float *g = NULL;
    float *gbasis;
    float *gx1, *gx2, *gy1, *gy2, *gz1, *gz2;
    float *k = NULL;
    float *kx, *ky;
    float kr;
    float *grew[3] = {NULL,NULL,NULL};
    float *g3 = NULL;
    float fovtemp;
    float kmax = 5/res;
    float gmaxRead;

    int Ndes = Tdes/dt;
    float kstart;
    int nstart;
    int n;
    int m;
    int Llo, Lhi, L;
    int nrew = 200;

    int glentemp;
    int klentemp;
    int ndim = 2;
    int ndimOut;
    int d;

    float phi;

    /* calculate maximum readout gradient */
    if(fovFilt)
	    gmaxRead = min(maxFilterGrad(fovFilt, dt), gmax);
    else
	    gmaxRead = gmax;

    Ndes += Ndes%2;
    if(gpre==NULL)
    {
		ndimGpre = ndim;
        ndimOut = ndim;
    }
    else
        ndimOut = ndimGpre;

        Llo = 1;

		Lhi = 0;
	   for(n=0; n<vd->nsteps; n++)
	   {
           /*kr = kmax*vd->kn[n];*/
           kr = vd->kr[n];
            getFovf(vd, &fovtemp);
		   if(sptype==spFERMAT)
			fovtemp = calcFovFermatFloret(fovtemp, kr, floretAngle);

			Lhi = max(Lhi, 2.0f*PI*kr*fovtemp);
	   }

       if(kr<kmax)
       {
           kr = kmax;
           getFovf(vd, &fovtemp);
          if(sptype==spFERMAT)
           fovtemp = calcFovFermatFloret(fovtemp, kr, floretAngle);

           Lhi = max(Lhi, 2.0f*PI*kr*fovtemp);

       }

        L = (Llo+Lhi)/2.0f;
        m = 0;

        if(!Tdes)
        {
            Llo = *nintl;
            Lhi = *nintl+2;
            L = *nintl;
            Ndes = MAX_GRAD_LENGTH;
        }

        printf("Number of desired points:\t%d\n", Ndes);

        gbasis = (float*)calloc(ndimOut*Ndes, sizeof(float));

        while((Lhi-Llo)>1 && m<50)
        {
            if(g)
            {
                free(g);
                g = NULL;
            }

            if(k)
            {
                free(k);
                k = NULL;
            }

		  if(genspiral(res, vd, L, Ndes, dt, sptype, floretAngle, gmaxRead, smax, &kstart, &nstart, &k, &g, &klentemp))
            {
                for(d=0; d<ndimGpre; d++)
                    if(grew[d])
                    {
                        free(grew[d]);
                        grew[d] = NULL;
                    }

                if(gpre)
                {
                    if(g3)
                    {
                        free(g3);
                        g3 = NULL;
                    }
                    n = nptsGpre+klentemp;
                    g3 = (float*)calloc(3*n, sizeof(float));
                    for(d=0; d<ndimGpre; d++)
                        memcpy(&g3[d*n], &gpre[d*nptsGpre], nptsGpre*sizeof(float));
                    for(d=0; d<ndim; d++)
                        memcpy(&g3[d*n+nptsGpre], &g[d*klentemp], klentemp*sizeof(float));
                    gx1 = g3;
                    gy1 = &g3[n];
                    if(ndimGpre>2)
                        gz1 = &g3[2*n];
                    else
                        gz1 = NULL;
				rewindTraj3(gx1, gy1, gz1, n, dt, 0, gmax, smax, &grew[0], &grew[1], &grew[2], &glentemp);
                }
                else
                {
				rewindTraj3(g, &g[klentemp], NULL, klentemp, dt, 0, gmax, smax, &grew[0], &grew[1], NULL, &glentemp);
                    nrew = glentemp - klentemp;
                }

                if(glentemp>Ndes)
                {
                    Llo = L;
                }
                else
                {
                    Lhi = L;
                    for(d=0; d<ndimOut; d++)
                        memcpy(&gbasis[d*Ndes], grew[d], glentemp*sizeof(float));

                    *Nk = klentemp;
                    *Ng = glentemp;
                    *nintl = L;
                }
            }
            else
                Llo = L;

            printf(" %d interleaves %d pts\n", L, glentemp);
            L = (Llo+Lhi)/2;

            if(!Tdes)
                Lhi = Llo;

            m++;
        }

    gx1 = gbasis;
    gy1 = &gbasis[Ndes];
    if(ndimOut>2)
        gz1 = &gbasis[2*Ndes];

    if(doOutputAll)
	    L = *nintl;
	else
		L = 1;

    *gAll = (float*)malloc(ndimOut**Ng*L*sizeof(float));
    *kAll = (float*)malloc(ndimOut**Nk*L*sizeof(float));
    if(dcfAll)
	    *dcfAll = (float*)malloc(*Nk*L*sizeof(float));

    if(gpre)
        m = nptsGpre + *Nk;
    else
        m = *Nk;

    for(n=0; n<L; n++)
    {
        phi = (2.0f*PI*n)/(*nintl);

        gx2 = &(*gAll)[ndimOut*n**Ng];
        gy2 = &(*gAll)[(ndimOut*n+1)**Ng];
        if(ndimOut>2)
            gz2 = &(*gAll)[(ndimOut*n+2)**Ng];

        kx = &(*kAll)[ndimOut*n**Nk];
        ky = &(*kAll)[(ndimOut*n+1)**Nk];

        memcpy(gx2, gx1, *Ng*sizeof(float));
        memcpy(gy2, gy1, *Ng*sizeof(float));
        scalecomplex(gx2, gy2, cos(phi), sin(phi), *Ng);

        if(gpre)
        {
            grad2k1(&gx2[nptsGpre], kx, dt, *Nk);
            grad2k1(&gy2[nptsGpre], ky, dt, *Nk);
        }
        else
        {
            grad2k1(gx2, kx, dt, *Nk);
            grad2k1(gy2, ky, dt, *Nk);
        }
        if(gpre && ndimOut>2)
        {
            memcpy(gz2, gz1, *Ng*sizeof(float));
            grad2k1(&gz2[nptsGpre], &(*kAll)[(ndimOut*n+2)**Nk], dt, *Nk);
        }

		if(dcfAll)
		  calcSpiralDcf(&gx2[nptsGpre], &gy2[nptsGpre], kx, ky, *Nk, &(*dcfAll)[n**Nk]);
    }

    return;
}

void makeSpiralTTraj(struct Vd *vd, float res, enum SpiralType sptype, float floretAngle, float Tdes, float dt, float fovFilt, float gmax, float smax, enum TrajStorage tst, struct Trajectory *spiral)
{
    float *g=NULL, *k=NULL, *dcf=NULL;
    int l;
    int d;

    int ng, nk;

    adjustRes(vd->fov[0], spiral->imgN, &res);

    makeSpiralT(vd, res, Tdes, dt, sptype, floretAngle, fovFilt, gmax, smax, &nk, &ng, &spiral->ninter, NULL, 0, 0, 1, &g, &k, &dcf);

    spiral->npts = ng;
    spiral->acqLength = nk;
    spiral->trajType = tjSPIRAL;
    spiral->dt = dt;
    getFov(vd, 0, spiral->FOV);
    spiral->FOV[1] = spiral->FOV[0];
    getFovf(vd, &spiral->FOV[2]);
    spiral->kmax = 5/res;
    spiral->max_G = gmax;
    spiral->max_slew = smax;
    spiral->res[0] = res;
    spiral->res[1] = res;
    spiral->imgN[1] = spiral->imgN[0];
    spiral->nbasis = 1;
    spiral->storage = tst;
    spiral->naxes = 2;

    if(tst==tstFULL)
    {
        spiral->waveformsG = matrix3f(spiral->naxes, spiral->ninter, spiral->npts);
        spiral->ks = matrix3f(spiral->naxes, spiral->ninter, spiral->acqLength);
        spiral->dcf = matrix2f(spiral->ninter, spiral->acqLength);

        for(l=0; l<spiral->ninter; l++)
        {
            for(d=0; d<spiral->naxes; d++)
            {
                memcpy(spiral->waveformsG[d][l], &g[(l*spiral->naxes+d)*spiral->npts], spiral->npts*sizeof(float));
                memcpy(spiral->ks[d][l], &k[(l*spiral->naxes+d)*spiral->acqLength], spiral->acqLength*sizeof(float));
            }
            memcpy(spiral->dcf[l], &dcf[l*spiral->acqLength], spiral->acqLength*sizeof(float));
        }
    }
    else
    {
        spiral->waveformsG = matrix3f(spiral->naxes, spiral->nbasis, spiral->npts);
        spiral->ks = matrix3f(spiral->naxes, spiral->nbasis, spiral->acqLength);
        spiral->dcf = matrix2f(spiral->nbasis, spiral->acqLength);;

        for(l=0; l<spiral->nbasis; l++)
        {
            for(d=0; d<spiral->naxes; d++)
            {
                memcpy(spiral->waveformsG[d][l], &g[(l*spiral->naxes+d)*spiral->npts], spiral->npts*sizeof(float));
                grad2k1(spiral->waveformsG[d][l], spiral->ks[d][l], spiral->dt, spiral->acqLength);
            }
            memcpy(spiral->dcf[l], &dcf[l*spiral->acqLength], spiral->acqLength*sizeof(float));
        }
    }

    if(g)
        free(g);
    if(k)
        free(k);
    if(dcf)
        free(dcf);

    return;
}

void makeArchSpiralTTraj(struct Vd *vd, float res, float Tdes, float dt, float fovFilt, float gmax, float smax, enum TrajStorage tst, struct Trajectory *spiral)
{
    makeSpiralTTraj(vd, res, spARCH, 0, Tdes, dt, fovFilt, gmax, smax, tst, spiral);

    return;
}

void makeArchSpiralTA(float fov, float res, float Tdes, float dt, float fovFilt, float gmax, float smax)
{
	struct Vd *vd;
	struct Trajectory spiral;


    makeSpiralTTraj(vd, res, spARCH, 0, Tdes, dt, fovFilt, gmax, smax, tstFULL, &spiral);

    return;
}

/*End of NOSPIRALGEN*/

