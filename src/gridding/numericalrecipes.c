#include "numericalrecipes.h"

#include <math.h>
#include <stdio.h>

/*
courtesy of NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
Copyright (C) 1988-1992 by Cambridge University Press
*/
float besseli(int n, float x)
{
	float ax,ans;
	double y;		/* Accumulate polynomials in double precision */

	switch(n)
	{
		case 0:
			if(x==0)
				ans = 1;
			else if ((ax=fabs(x)) < 3.75)  /* Polynomial t. */
			{
				y = x/3.75;
				y *= y;
				ans = 1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
			}
			else
			{
				y = 3.75/ax;
				ans = (exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))))))));
			}
			break;
		default:
			printf("Invalid bessel order entered\n");
	}
	return ans;
}
