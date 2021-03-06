#include "phantom.h"

#include "shape.h"
extern "C"
{
#include "arrayops.h"
}

#include <algorithm>
/*
 * SheppLogan3D.java
 *
 * Created on January 7, 2007, 8:46 PM
 */

/**
 *
 * Three-dimensional Shepp-Logan Phantom in both the Fourier and image domains.
 *
 * This is a class called SheppLogan3D. It can be used to generate Fourier domain
 * signal (or k-space) as well as image domain signal based on a 3-dimensional
 * analytical phantom in both domains. Please refer to the source code or
 * the article referenced below for further information.
 *
 * <br>
 * <br> The image below is a 3D rendering of the phantom in the image domain</p>
 * <p><p align="center"><img src="threed.png" width="285" height="345"></p><br></p>
 *
 *
 * <br>
 * The reconstructed images from the Fourier domain signals are shown below as animated gif:
 *
 * <p><p align="center"><img src="reconmovie.gif" width="345" height="345"></p><br></p>
 *
 * <br>
 *
 *
 * <br>
 * Please refer to
 * <br>
 * Koay CG, Sarlls JE, &#214zarslan E.
 * Three Dimensional Analytical Magnetic Resonance Imaging Phantom in the Fourier Domain. Magn Reson Med. 58: 430-436 (2007)
 * <br>
 * for further information.
 *
 * <br>
 * @see <a href=http://dx.doi.org/10.1002/mrm.21292>Ref</a>
 * @author  Cheng Guan Koay
 * @since 07/25/2007
 *
 * @version &#8722&#8734.
 */

/** Creates a new instance of SheppLogan3D */
Phantom::Phantom(std::vector<float> fieldOfView)
{
	float scale = 0;
	for(size_t d=0; d<fieldOfView.size(); d++)
		scale = std::max(fieldOfView[d], scale);
	scale *= .5;

	if(fieldOfView.size()==3)
	{
		m_shapes.push_back(Shape(Shape::Ellipsoid, 0,       0,       0,     0.69*scale,    0.92*scale,     0.9*scale,              0,      0,    0,      2.));
    m_shapes.push_back(Shape(Shape::Ellipsoid, 0,       0,       0,   0.6624*scale,   0.874*scale,    0.88*scale,              0,      0,    0,    -0.8f));
    m_shapes.push_back(Shape(Shape::Ellipsoid, -0.22*scale,      0.,   -0.25*scale,     0.41*scale,    0.16*scale,    0.21*scale, static_cast<float>(3*M_PI)/5,      0,    0,    -0.2f));
    m_shapes.push_back(Shape(Shape::Ellipsoid, 0.22*scale,      0.,   -0.25*scale,     0.31*scale,    0.11*scale,    0.22*scale, static_cast<float>(2*M_PI)/5,      0,    0,    -0.2f));
    m_shapes.push_back(Shape(Shape::Ellipsoid, 0,    0.35*scale,   -0.25*scale,     0.21*scale,    0.25*scale,     0.5*scale,              0,      0,    0,     0.2f));
    m_shapes.push_back(Shape(Shape::Ellipsoid, 0,     0.1*scale,   -0.25*scale,    0.046*scale,   0.046*scale,   0.046*scale,              0,      0,    0,     0.2f));
    m_shapes.push_back(Shape(Shape::Ellipsoid, -0.08*scale,   -0.65*scale,   -0.25*scale,    0.046*scale,   0.023*scale,    0.02*scale,              0,      0,    0,     0.1f));
    m_shapes.push_back(Shape(Shape::Ellipsoid, 0.06*scale,   -0.65*scale,   -0.25*scale,    0.046*scale,   0.023*scale,    0.02*scale,              0,      0,    0,     0.1f));
    m_shapes.push_back(Shape(Shape::Ellipsoid, 0.06*scale,  -0.105*scale,   0.625*scale,    0.056*scale,    0.04*scale,     0.1*scale,     static_cast<float>(M_PI/2),      0,    0,     0.2f));
    m_shapes.push_back(Shape(Shape::Ellipsoid, 0.,     0.1*scale,   0.625*scale,    0.056*scale,   0.056*scale,     0.1*scale,     static_cast<float>(M_PI/2),      0,    0,    -0.2f));
	}
	else
	{
    m_shapes.push_back(Shape(Shape::Ellipse, 0,       0,        0.92*scale,    0.69*scale,     static_cast<float>(M_PI_2),       2));
    m_shapes.push_back(Shape(Shape::Ellipse,        0, -0.0184*scale,        0.874*scale,   0.6624*scale,   static_cast<float>(M_PI_2),     -0.8f));
    m_shapes.push_back(Shape(Shape::Ellipse,     0.22*scale,      0.,        0.31f,    0.11*scale,  static_cast<float>(2*M_PI/5),     -0.2f));
    m_shapes.push_back(Shape(Shape::Ellipse,    -0.22*scale,      0.,        0.41*scale,    0.16f,  static_cast<float>(3*M_PI/5),     -0.2f));
    m_shapes.push_back(Shape(Shape::Ellipse,        0,    0.35*scale,        0.25*scale,    0.21*scale,     static_cast<float>(M_PI_2),     -0.2f));
    m_shapes.push_back(Shape(Shape::Ellipse,        0,     0.1*scale,        0.046*scale,   0.046*scale,              0,      0.1f));
    m_shapes.push_back(Shape(Shape::Ellipse,       0.,    -0.1*scale,        0.046*scale,   0.046*scale,              0,      0.1f));
    m_shapes.push_back(Shape(Shape::Ellipse,    -0.08*scale,   -0.605*scale,       0.046*scale,   0.023*scale,              0,      0.1f));
    m_shapes.push_back(Shape(Shape::Ellipse,       0.,   -0.605*scale,       0.023*scale,   0.023*scale,              0,      0.1f));
    m_shapes.push_back(Shape(Shape::Ellipse,     0.06*scale,   -0.605*scale,       0.046*scale,   0.023*scale,    static_cast<float>(M_PI_2),      0.1f));

	}
}

/**
 *  User may add new ellipsoids and change their properties with this constructor.
 *
 *  @param ellipsoids is a set of ellipsoids describing the phantom
 *
 *  Please to the paper mentioned above for further information on the notations.
 *
 */
Phantom::Phantom(std::vector<Shape> shapes)
{
	for(size_t n=0; n<m_shapes.size(); n++)
		m_shapes.push_back(shapes[n]);
}


/**
  *  Given a list of position vectors, i.e. {{x1,y1,z1},{x2,y2,z2},...},
  *  the image domain signals at those locations are returned.
  *
  */
std::vector<float> Phantom::imageDomainSignal(const std::vector<float>& coordinates)
{
	 int points = coordinates.size()/3;
	 std::vector<float> data(points);

	 for(int i=0; i<points; i++)
	{
		 data[i] = imageDomainSignal(coordinates[3*i], coordinates[3*i+1], coordinates[3*i+2]);

	 }

	 return data;
}


/**
 * returning real value of the image intensity at (x,y,z).
 *
 */
float Phantom::imageDomainSignal(float x, float y, float z)
{

	float position[3] = {x,y,z};
	double signal = 0.0;

	for(size_t i=0; i<m_shapes.size(); i++){ // loop through each of the ellipsoids
		Shape& shape = m_shapes.at(i);
		float relativePosition[3];
		shape.relativePosition(position, relativePosition);
		float sum = 0.0;
		 for(int d=0; d<shape.dimensions(); d++)
		 {
			 float projection = relativePosition[d]/shape.principalAxis(d);
			sum += projection*projection;
		 }
		 signal += (sum<=1.0) ? shape.intensity() : 0;
	}

	return signal;
}

float Phantom::imageDomainSignal(float x, float y)
{
	return imageDomainSignal(x,y,0);
}
 /**
  *  Given a list of (kx,ky,kz), the k-space signals at those locations are returned.
  *  The return array is of dimension kList.length by 2.
  *  The first column of the array is the real part of the complex signal and the second is
  *  the imaginary part of the complex signal.
  */
 std::vector<complexFloat> Phantom::fourierDomainSignal(const std::vector<float>& coordinates)
 {
	 int points = coordinates.size();
	 std::vector<complexFloat> data(points);

	 for(int i=0; i<points; i++)
	 {
		 data[i] = fourierDomainSignal(coordinates[3*i], coordinates[3*i+1], coordinates[3*i+2]);
	 }

	 return data;
 }


/**
 * returning the complex signal evaluated at ( kx, ky, kz) in an array of length 2, i.e. {Re, Im}.
 */
complexFloat Phantom::fourierDomainSignal(const float& kx, const float& ky, const float& kz)
{
	float k[3] = {kx,ky,kz};
	float kRotated[3];
	float kProjected[3];

	complexFloat signal; // {Re, Im} , real and imaginary signals

	double arg = 0.0;
	for(size_t i=0; i<m_shapes.size(); i++)
	{
		Shape& shape = m_shapes.at(i);
		float principalAxesProduct = norm2(shape.principalAxes().data(), 3);

		shape.rotatedCoordinates(k, kRotated);
		multiplyfloats(shape.principalAxes().data(), kRotated, kProjected, 3);
		float K = norm2(kProjected, 3);

		 arg = 2.0 * M_PI * K;

		 if(K==0.0)
		 { // if K = 0
			 if( norm2(shape.displacement().data(),3)==0.0 )
			 { // if displacement vector is zero
         signal += (4.0f/3.0f)*M_PI* shape.intensity()*principalAxesProduct;
			 }
			 else
			 { // displacement vector is not zero
				 double kd = dot(k, shape.displacement().data(), 3);
				 double temp = (4./3.)*M_PI* shape.intensity()*principalAxesProduct;
				 // signal.real() += temp * cosf(2.0 * M_PI * kd);
				 // signal.imag() -= temp * sinf(2.0 * M_PI * kd);
                 signal += complexFloat(temp * std::cos(2 * M_PI * kd), -temp * std::sin(2 * M_PI * kd));
			 }

		 }else if (K<=0.002){  // if K<=0.002
			 if( norm2(shape.displacement().data(),3)==0.0 ){ // if displacement vector is zero

				 double temp = 4.1887902047863905 - 16.5366808961599*powf(K,2) + 23.315785507450016*powf(K,4);
         signal += shape.intensity()*principalAxesProduct*temp;

			 }else{  // if displacement vector is not zero
				 double kd = dot(k, shape.displacement().data(), 3);
				 double temp1 = 4.1887902047863905 - 16.5366808961599*powf(K,2) + 23.315785507450016*powf(K,4);
				 double temp2 = shape.intensity()*principalAxesProduct*temp1;

//				 signal.real() += temp2 * cosf(2.0 * M_PI * kd);
//				 signal.imag() -= temp2 * sinf(2.0 * M_PI * kd);
        signal += complexFloat(temp2 * std::cos(2 * M_PI * kd), -temp2 * std::sin(2 * M_PI * kd));
			 }
		 }
		 else
		 { // K>0.002
			 if( norm2(shape.displacement().data(),3)==0.0 ){ // if displacement vector is zero

         float temp = sinf(arg)-arg*cosf(arg);
            temp /= (2*std::pow(static_cast<float>(M_PI), 2) * std::pow(K, 3));

         signal += shape.intensity()*principalAxesProduct*temp;

			 }else{  // displacement vector is not zero
				 double kd = dot(k, shape.displacement().data(), 3);
         float temp = sinf(arg)-arg*cosf(arg);
            temp /= (2.0*powf(static_cast<float>(M_PI), 2)*powf(K,3));

						temp *= shape.intensity()*principalAxesProduct;

//				 signal.real() += temp * cosf(2.0 * M_PI * kd);
//				 signal.imag() -= temp * sinf(2.0 * M_PI * kd);
        signal += complexFloat(temp * std::cos(2.0 * M_PI * kd), -temp * std::sin(2.0 * M_PI * kd));
			 }
		 }
	}
  
	return signal;
}

complexFloat Phantom::fourierDomainSignal(float kx, float ky)
{
	float k[2] = {kx,ky};
	float kRotated[2];
	float kProjected[2];

	complexFloat signal; // {Re, Im} , real and imaginary signals

	double K = 0.0;
	double K2 = 0.0;
	double K4 = 0.0;

	for(size_t i=0; i<m_shapes.size(); i++)
	{
		Shape& shape = m_shapes.at(i);
		float principalAxesProduct = shape.principalAxis(0)*shape.principalAxis(1);

//		 K = norm2( MAT.multiply(dot(RT[i],k), ab[i]) );
		shape.rotatedCoordinates(k, kRotated);
		multiplyfloats(shape.principalAxes().data(), kRotated, kProjected, 2);
		K = norm2(kProjected, 2);

		 if(K==0.0){ // if K = 0
			 signal += M_PI* shape.intensity()*principalAxesProduct;

		 }else if (K<=0.002){  // if K<=0.002
			 if(norm2(shape.displacement().data(),2)==0.0)
			 {
					K2 = K*K;
					K4 = K2*K2;

					double BJ = 1. - 4.934802200544679*K2 + 8.117424252833533*K4;

					signal += M_PI* shape.intensity()*principalAxesProduct*BJ;

			 }else{
					K2 = K*K;
					K4 = K2*K2;

          const float BJ = 1. - 4.934802200544679*K2 + 8.117424252833533*K4;

          const float kd = dot(k, shape.displacement().data(), 2);
          const float arg2 = 2.0 * M_PI * kd;

          const float temp = M_PI* shape.intensity()*principalAxesProduct;
//					signal.real() += temp * std::cos(arg2)*BJ;
//					signal.imag() -= temp * std::sin(arg2)*BJ;
          signal += complexFloat(temp * std::cos(arg2)*BJ, -temp * std::sin(arg2)*BJ);
			 }
		 }
		 else
		 { // K>0.002
 if(norm2(shape.displacement().data(),2)==0.0)
			 {
				double arg = M_PI * K;
				double BJ = j1(2.0*arg)/arg;
        signal += M_PI * shape.intensity()*principalAxesProduct*BJ;
			 }
			 else
			 {
        const float arg = M_PI * K;
        const float BJ = j1(2.0*arg)/arg;
        const float kd = dot(k, shape.displacement().data(), 2);
        const float arg2 = 2.0 * M_PI * kd;

        const float temp = M_PI* shape.intensity()*principalAxesProduct;
//				signal.real() += temp * std::cos(arg2)*BJ;
//				signal.imag() -= temp * std::sin(arg2)*BJ;
        signal += complexFloat(temp * std::cos(arg2)*BJ, -temp * std::sin(arg2)*BJ);
			 }
		 }
	}

	return signal;
}
