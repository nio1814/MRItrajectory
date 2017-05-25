#include "ellipsoid.h"

extern "C"
{
#include "mathops.h"
}

#include <math.h>

void Rx(float t, float* rotationMatrix)
{
	rotationMatrix[0] = 1;
	rotationMatrix[1] = 0;
	rotationMatrix[2] = 0;
	rotationMatrix[3] = 0;
	rotationMatrix[4] = cosf(t);
	rotationMatrix[5] = -sinf(t);
	rotationMatrix[6] = 0;
	rotationMatrix[7] = sinf(t);
	rotationMatrix[8] = cosf(t);
}

void Ry(float t, float* rotationMatrix)
{
   rotationMatrix[0] = cosf(t);
   rotationMatrix[1] = 0;
   rotationMatrix[2] = sinf(t);
   rotationMatrix[3] = 0;
   rotationMatrix[4] = 1;
   rotationMatrix[5] = 0;
   rotationMatrix[6] = -sinf(t);
   rotationMatrix[7] = 0;
   rotationMatrix[8] = cosf(t);
}

void Rz(float t, float* rotationMatrix)
{
   rotationMatrix[0] = cosf(t);
   rotationMatrix[1] = -sinf(t);
   rotationMatrix[2] = 0;
   rotationMatrix[3] = sinf(t);
   rotationMatrix[4] = cosf(t);
   rotationMatrix[5] = 0;
   rotationMatrix[6] = 0;
   rotationMatrix[7] = 0;
   rotationMatrix[8] = 1;
}

Ellipsoid::Ellipsoid(float x, float y, float z, float xAxisLength, float yAxisLength, float zAxisLength, float phi, float theta, float psi, float intensity)
{
	m_displacement[0] = x;
	m_displacement[1] = y;
	m_displacement[2] = z;

	m_principalAxes[0] = xAxisLength;
	m_principalAxes[1] = yAxisLength;
	m_principalAxes[2] = zAxisLength;

	float rotationMatrixX[9];
	Rx(theta, rotationMatrixX);
	float rotationMatrixY[9];
	Ry(psi, rotationMatrixY);
	float rotationMatrixZ[9];
	Rz(phi, rotationMatrixZ);

	float Rzy[9];
	matrixMultiply(rotationMatrixZ, 3, 3, rotationMatrixY, 3, Rzy);

	matrixMultiply(Rzy, 3, 3, rotationMatrixX, 3, m_rotationMatrix);

	m_intensity = intensity;
}

void Ellipsoid::relativePosition(const float* position, float* positionRelative)
{
	for(int d=0; d<3; d++)
	{
		float offset = position[d]-m_displacement[d];
		positionRelative[d] = 0;
		for(int n=0; n<3; n++)
			positionRelative[d] += m_rotationMatrix[3*n+d]*offset;
	}
}

void Ellipsoid::rotatedCoordinates(const float *coordinates, float* coordinatesRotated)
{
	for(int m=0; m<3; m++)
	{
		coordinatesRotated[m] = 0;
		for(int n=0; n<3; n++)
			coordinatesRotated[m] += m_rotationMatrix[3*m+n]*coordinates[m];
	}
}
