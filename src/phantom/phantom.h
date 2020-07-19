#ifndef PHANTOM_H
#define PHANTOM_H

#include "shape.h"

#include <complex>

using complexFloat = std::complex<float>;

class Phantom
{
public:
	Phantom(std::vector<float> fieldOfView);
	Phantom(std::vector<Shape> shapes);

	std::vector<float> imageDomainSignal(const std::vector<float>& coordinates);
	float imageDomainSignal(float x, float y, float z);
	float imageDomainSignal(float x, float y);

	std::vector<complexFloat> fourierDomainSignal(const std::vector<float>& coordinates);

  /*!
   * \brief Calculate the complex fourier phantom signal at a k-space location.
   * \param kx
   * \param ky
   * \param kz
   */
  complexFloat fourierDomainSignal(const float& kx, const float& ky, const float& kz);
	complexFloat fourierDomainSignal(float kx, float ky);

private:

	std::vector<Shape> m_shapes;
};

#endif // PHANTOM_H
