#ifndef PHANTOM_H
#define PHANTOM_H

#include "shape.h"

#include <complex>

typedef std::complex<float> complexFloat;

class Phantom
{
public:
	Phantom(std::vector<float> fieldOfView);
	Phantom(std::vector<Shape> shapes);

	std::vector<float> imageDomainSignal(const std::vector<float>& coordinates);
	float imageDomainSignal(float x, float y, float z);
	float imageDomainSignal(float x, float y);
private:

	std::vector<complexFloat> fourierDomainSignal(const std::vector<float>& coordinates);
	complexFloat fourierDomainSignal(float kx, float ky, float kz);

	std::vector<Shape> m_shapes;
};

#endif // PHANTOM_H
