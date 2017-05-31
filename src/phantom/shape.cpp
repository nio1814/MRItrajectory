#include "shape.h"

extern "C"
{
#include "mathops.h"
}

#include <cmath>

void Rx(float t, float* rotationMatrix)
{
	rotationMatrix[0] = 1;
	rotationMatrix[1] = 0;
	rotationMatrix[2] = 0;
	rotationMatrix[3] = 0;
	rotationMatrix[4] = std::cos(t);
	rotationMatrix[5] = -sinf(t);
	rotationMatrix[6] = 0;
	rotationMatrix[7] = sinf(t);
	rotationMatrix[8] = std::cos(t);
}

void Ry(float t, float* rotationMatrix)
{
   rotationMatrix[0] = std::cos(t);
   rotationMatrix[1] = 0;
   rotationMatrix[2] = sinf(t);
   rotationMatrix[3] = 0;
   rotationMatrix[4] = 1;
   rotationMatrix[5] = 0;
   rotationMatrix[6] = -sinf(t);
   rotationMatrix[7] = 0;
   rotationMatrix[8] = std::cos(t);
}

void Rz(float t, float* rotationMatrix)
{
   rotationMatrix[0] = std::cos(t);
   rotationMatrix[1] = -sinf(t);
   rotationMatrix[2] = 0;
   rotationMatrix[3] = sinf(t);
   rotationMatrix[4] = std::cos(t);
   rotationMatrix[5] = 0;
   rotationMatrix[6] = 0;
   rotationMatrix[7] = 0;
   rotationMatrix[8] = 1;
}

Shape::Shape(Type type, float x, float y, float xAxisLength, float yAxisLength, float theta, float intensity) :
	m_intensity(intensity), m_type(type)
{
	m_displacement.resize(2);
	m_displacement[0] = x;
	m_displacement[1] = y;

	m_principalAxes.resize(2);
	m_principalAxes[0] = xAxisLength;
	m_principalAxes[1] = yAxisLength;

	m_rotationMatrix.resize(4);
	m_rotationMatrix[0] = std::cos(theta);
	m_rotationMatrix[1] = -std::sin(theta);
	m_rotationMatrix[2] = std::sin(theta);
	m_rotationMatrix[3] = std::cos(theta);
}

Shape::Shape(Type type, float x, float y, float z, float xAxisLength, float yAxisLength, float zAxisLength, float phi, float theta, float psi, float intensity) : m_type(type)
{
	m_displacement.resize(3);
	m_displacement[0] = x;
	m_displacement[1] = y;
	m_displacement[2] = z;

	m_principalAxes.resize(3);
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

	m_rotationMatrix.resize(9);
	matrixMultiply(Rzy, 3, 3, rotationMatrixX, 3, m_rotationMatrix.data());

	m_intensity = intensity;
}


float Shape::principalAxis(int axis)
{
	return m_principalAxes[axis];
}

std::vector<float> Shape::principalAxes()
{
	return m_principalAxes;
}

std::vector<float> Shape::displacement()
{
	return m_displacement;
}

float Shape::intensity()
{
	return m_intensity;
}

int Shape::dimensions()
{
	return m_principalAxes.size();
}

void Shape::relativePosition(const float* position, float* positionRelative)
{
	for(int d=0; d<dimensions(); d++)
	{
		float offset = position[d]-m_displacement[d];
		positionRelative[d] = 0;
		for(int n=0; n<dimensions(); n++)
			positionRelative[d] += m_rotationMatrix[dimensions()*n+d]*offset;
	}
}

void Shape::rotatedCoordinates(const float *coordinates, float* coordinatesRotated)
{
	for(int m=0; m<dimensions(); m++)
	{
		coordinatesRotated[m] = 0;
		for(int n=0; n<dimensions(); n++)
			coordinatesRotated[m] += m_rotationMatrix[dimensions()*m+n]*coordinates[m];
	}
}
