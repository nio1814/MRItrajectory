#ifndef VD_H
#define VD_H

enum VariableDensityFunction{VariableDensityPolynomial};

#define MAX_VARIABLE_DENSITY_STEPS 10

struct VariableDensityStep
{
	float kr;
	float param;
    enum VariableDensityFunction function;
    float scale;
};

struct VariableDensity
{
	int steps;
	struct VariableDensityStep step[MAX_VARIABLE_DENSITY_STEPS];
	float fovcomp;
	float *kcomp;
	int nkcomp;
	int compSelect;
};

void initializeVariableDensity(struct VariableDensity *v);

void copyVariableDensity(const struct VariableDensity *from, struct VariableDensity *to);

//void deleteVd(struct Vd *v);
float getScale(const struct VariableDensity *v, float kr);
float getFinalScale(const struct VariableDensity *v);
//void getFovComp(struct Vd *v, float kr, float *fov);

void addVariableDensityStep(struct VariableDensity *v, enum VariableDensityFunction function, float kr, float functionParameter, float scale);
int writeVariableDensity(const char* filename, const struct VariableDensity *v, float krMax, int points);
//void copyVdArrs(struct Vd *vdfrom, struct Vd *vdto);

//void calcAngleComp(struct Vd *v, float fov, float res);

///**
//  \brief	Anisoptropic field of view
//  */
//void calcAngleCompa(struct Vd *v, float fovr, float fovp, float kmax);

//void calcAngleCompb(struct Vd *v, float td, float fovx, float fovy, float res, int isConstAng);

void getFieldOfView(const struct VariableDensity* variableDensity, float kr, const float *fieldOfViewInitial, float *fieldOfViewKr, int dimensions);

void getFinalFieldOfView(const struct VariableDensity *variableDensity, const float *initialFieldOfView, float *finalFieldOfView, int dimensions);

#endif

