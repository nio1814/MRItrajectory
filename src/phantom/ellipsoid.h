#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include <vector>

struct Ellipsoid
{
	Ellipsoid(float x, float y, float z, float xAxisLength, float yAxisLength, float zAxisLength, float phi, float theta, float psi, float intensity);

	void relativePosition(const float *position, float* positionRelative);

	void rotatedCoordinates(const float *coordinates, float* coordinatesRotated);

	float m_rotationMatrix[9];	//rotation matrices. T denotes matrix transpose.
	float m_displacement[3]; //displacement vectors.
	float m_principalAxes[3]; //the length of the principal axes.
	float m_intensity; // signal intensity.
};

#endif // ELLIPSOID_H
