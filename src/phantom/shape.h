#ifndef SHAPE_H
#define SHAPE_H

#include <vector>

class Shape
{
public:
	enum Type{Ellipsoid, Ellipse};

	Shape(Type type, float x, float y, float xAxisLength, float yAxisLength, float theta, float intensity);
	Shape(Type type, float x, float y, float z, float xAxisLength, float yAxisLength, float zAxisLength, float phi, float theta, float psi, float intensity);

	void relativePosition(const float *position, float* positionRelative);

	void rotatedCoordinates(const float *coordinates, float* coordinatesRotated);

	float principalAxis(int axis);
	std::vector<float> principalAxes();
	std::vector<float> displacement();
	float intensity();

	int dimensions();

private:
	std::vector<float> m_rotationMatrix;	//rotation matrices. T denotes matrix transpose.
	std::vector<float>  m_displacement; //displacement vectors.
	std::vector<float>  m_principalAxes; //the length of the principal axes.
	float  m_intensity; // signal intensity.
	Type m_type;
};

#endif // SHAPE_H
