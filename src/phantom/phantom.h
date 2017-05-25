#ifndef PHANTOM_H
#define PHANTOM_H

#include "ellipsoid.h"

#include <vector>
#include <complex>

typedef std::complex<float> complexFloat;

class Phantom
{
public:
	Phantom(std::vector<float> fieldOfView);
	Phantom(std::vector<Ellipsoid> ellipsoids);

	std::vector<float> imageDomainSignal(const std::vector<float>& coordinates);
	float imageDomainSignal(float x, float y, float z);
private:

	std::vector<complexFloat> fourierDomainSignal(const std::vector<float>& coordinates);
	complexFloat fourierDomainSignal(float kx, float ky, float kz);

	std::vector<Ellipsoid> m_ellipsoids;
};

#endif // PHANTOM_H
