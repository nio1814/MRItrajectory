#ifndef SPIRAL_H
#define SPIRAL_H

#include "trajectory_s.h"
#include "vd.h"

#define SPIRALDT 1e-6
#define MAX_GRAD_LENGTH	8192

/**
 \ingroup	trajectory
 \defgroup	spiral	Spiral
 @{ 
*/

/**
 \defgroup	spiral	Spirals
 @{
*/

enum SpiralDirection{spOUT, spINOUT};
enum spMETHOD{spmVDS, spmMINTIME};
enum SpiralType{spARCH, spFERMAT};

void printSpiral2Info(struct Trajectory* spiral);

/**
 Make instruction waveforms based on first interleaf
 */
void makeInterleavesI(struct Trajectory* spiral, int doWrite);

/**
 Make a set of interleaves
 \param[in]	wavesInit	Basis waveforms
 \param[in]	ninter	Number of interleaves
 \param[in]	npts	Number of points in waveform
 \param[out]	wavesAll	Interleaved waveforms
 */
void makeInterleaves(float** wavesInit, int ninter, int npts, float*** wavesAll);

/**
 \param[in]	emodek	Endian mode of k-space trajectory file
*/
int spiralLoadTrajRth(char* ksfile, char* wavfile, struct Trajectory *traj, char emodek, char emodew);

/**
 \param[in]	doK	Make k-space coordinate interleaves
 \param[in]	doG	Make gradient waveform interleaves
 \param[in]	doI	Make instruction amplitude interleaves
*/
void makeSpiraInterleavesTraj(struct Trajectory* spiral, int doK, int doG, int doI);

void kspaceReconFile(struct Trajectory* spiral);

float calcc1Spiral(float fov, float nintl);
float calcGtwistSpiral(float kr, float fr, float nintl);

/**
 \brief	Calculate the sampling density compensation weights based on curvature and spacing of samples
 \param[in]	gx	x gradient values (G/cm)
 \param[in]	gy	y gradient values (G/cm)
*/
void calcSpiralDcf(float *gx, float *gy, float *kx, float *ky, int rolen, float *denscomp);


/**
 \brief	Make spiral trajectory based on Jim Pipe's code
*/
/*void makeSpiralJgpTraj(float fovi, float fovf, float res, int nintl, float smax, float gmax, float sysgmax, int Ts_us, struct Trajectory *spiral);*/
void makeSpiralJgpTraj(float fovi, float fovf, float vdStart, float vdEnd, int vdType, float vdParam, float res, float Tdes, int nintl, float smax, float gmax, float sysgmax, int Ts_us, struct Trajectory *spiral);

void makeSpiral(struct Vd *vd, float res, float Tdes, float dt, float gmax, float smax, int *Nk, int *Ng, int *nintl, float **gAll, float **kAll, float **dcfAll);

/**
 \param[in]	res	Resolution (mm)
 \param[in]	Tdes	Waveform duration (sec)
 \param[in]	dt	Sampling period (sec)
 \param[in]	gmax	Gradient amplitude limit (G/cm)
 \param[in]	smax	Slew rate limit (G/cm/s)
*/
void makeSpiralT(struct Vd *vd, float res, float Tdes, float dt, enum SpiralType sptype, float floretAngle, float fovfilt, float gmax, float smax, int *Nk, int *Ng, int *nintl, float *gpre, int nptsGpre, int ndimGpre, int doOutputAll, float **gAll, float **kAll, float **dcfAll);
void makeArchSpiralTTraj(struct Vd *vd, float res, float Tdes, float dt, float fovFilt, float gmax, float smax, enum TrajStorage tst, struct Trajectory *spiral);
void makeSpiralTTraj(struct Vd *vd, float res, enum SpiralType sptype, float floretAngle, float Tdes, float dt, float fovFilt, float gmax, float smax, enum TrajStorage tst, struct Trajectory *spiral);

/** @} */

#endif

