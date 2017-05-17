#include "variabledensity.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void initializeVariableDensity(struct VariableDensity *v)
{
	v->steps = 0;
	addVariableDensityStep(v, VariableDensityPolynomial, 0, 0, 1);

	v->kcomp = NULL;
	v->compSelect = -1;

	return;
}

void copyVariableDensity(const struct VariableDensity *from, struct VariableDensity *to)
{
	int s;
	int totalSteps = from->steps;

	for(s=0; s<totalSteps; s++)
		to->step[s] = from->step[s];
	to->steps = from->steps;
	to->compSelect = from->compSelect;
}

void addVariableDensityStep(struct VariableDensity *v, enum VariableDensityFunction function, float kr, float functionParameter, float scale)
{
	int s=0;

    while(s<v->steps && (kr > v->step[s].kr))
		s++;

	if(s==v->steps)
	{
        v->step[s].kr = kr;
		v->step[s].param = functionParameter;
		v->step[s].function = function;
        /*memcpy(&v->fov[v->ndim*s], fov, v->ndim*sizeof(float));*/
        v->step[s].scale = scale;
		v->steps++;
	}
}

float getFinalScale(const struct VariableDensity *v)
{
//	memcpy(finalFieldOfView, v->fieldOfView, v->dimensions*sizeof(float));
	return  v->step[v->steps-1].scale;
}


float getScale(const struct VariableDensity *v, float kr)
{
	int s=0;
	float deltaKr;
	float deltaScale;
	float scale;

	while((s<v->steps-1) && (kr >= v->step[s].kr))
		s++;

	if((s==v->steps-1) && (kr >= v->step[s].kr))
	{
		/*for(n=0; n<v->ndim; n++)
			fieldOfView[n] = v->scale[s]*v->fieldOfView[n];*/
		  scale = v->step[s].scale;
	}
	else
	{
		s--;
		deltaKr = v->step[s+1].kr - v->step[s].kr;
		kr -= v->step[s].kr;
		deltaScale = v->step[s+1].scale - v->step[s].scale;

		switch(v->step[s].function)
		{
			case VariableDensityPolynomial:
				if(v->step[s].param>0)
					scale = v->step[s].scale + deltaScale*pow(kr/deltaKr, v->step[s].param);
				else if(v->step[s].param<0)
					scale = v->step[s+1].scale - deltaScale*pow(1.0f-kr/deltaKr, -v->step[s].param);
				else
					scale = 1;
				break;
			default:
				fprintf(stderr, "Error getScale: Invalid vd function %d\n", v->step[s].function);
				break;
		}
	}

//	if(v->compSelect>-1)
//		fieldOfViewScale *= getfieldOfViewCompScale(v, kSpaceRadius);
	return scale;
}

void printVariableDensity(const struct VariableDensity *v)
{
	int s;
	for(s=0; s<v->steps; s++)
		printf("%f\t%f\n", v->step[s].kr, v->step[s].scale);
}

int writeVariableDensity(const char* filename, const struct VariableDensity *v, float krMax, int points)
{
	FILE* file;
	int n;
	float kr;
	float scale;

	file = fopen(filename, "wb");
	if(!file)
	{
	   fprintf(stderr, "Error opening %s for read", filename);
	   return 1;
	}

	for(n=0; n<points; n++)
	{
		kr = n*krMax/points;
		fwrite(&kr, sizeof(float), 1, file);
		scale = getScale(v, kr);
		fwrite(&scale, sizeof(float), 1, file);
	}
	fclose(file);

	return 0;
}
