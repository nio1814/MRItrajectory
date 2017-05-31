#include "phantomtest.h"

#include "phantom.h"

#include <QtTest/QTest>
#include <QtMath>
#include <QTextStream>

void writeOctaveFileFloat(QString filename, const std::vector<float>& data, std::vector<int> imageSize)
{
	QFile file(filename);
	if(!file.open(QIODevice::WriteOnly))
	{
		qWarning("Failed to open %s", filename.toStdString().c_str());
		return;
	}

	QTextStream stream(&file);
	stream << "# Phantom Test\n";
	stream << "# name: m\n";
	stream << "# type: matrix\n";
	stream << "# ndims: " << imageSize.size() << endl;
	for(size_t d=0; d<imageSize.size(); d++)
		stream << " " << imageSize[d];
	endl(stream);

	for (size_t n=0; n<data.size(); n++)
		stream << " " << data[n] << endl;

	file.close();
}

void writeOctaveFileFloat(QString filename, const std::vector<complexFloat>& data, std::vector<int> imageSize)
{
	QFile file(filename);
	if(!file.open(QIODevice::WriteOnly))
	{
		qWarning("Failed to open %s", filename.toStdString().c_str());
		return;
	}

	QTextStream stream(&file);
	stream << "# Phantom Test\n";
	stream << "# name: m	\n";
	stream << "# type: complex matrix\n";
	stream << "# ndims: " << imageSize.size() << endl;
	for(size_t d=0; d<imageSize.size(); d++)
		stream << " " << imageSize[d];
	endl(stream);

	for (size_t n=0; n<data.size(); n++)
		stream << " (" << data[n].real() << ',' << data[n].imag() << ")\n";

	file.close();
}


void PhantomTest::testImage_data()
{
	QTest::addColumn<std::vector<float> >("fieldOfView");
	QTest::addColumn<std::vector<float> >("spatialResolution");

	std::vector<float> fieldOfView(3);
	std::vector<float> spatialResolution(3);
	for(int d=0; d<3; d++)
	{
		fieldOfView[d] = 28;
		spatialResolution[d] = 4;
	}
//	QTest::newRow("3D") << fieldOfView << spatialResolution;

	fieldOfView.resize(2);
	for(int d=0; d<2; d++)
		spatialResolution[d] = 1;

	QTest::newRow("2D") << fieldOfView << spatialResolution;
}

void PhantomTest::testImage()
{
	QFETCH(std::vector<float>, fieldOfView);
	QFETCH(std::vector<float>, spatialResolution);

	int dimensions = fieldOfView.size();
	std::vector<int> imageSize(dimensions);
	int points = 1;
	for(int d=0; d<dimensions; d++)
	{
		imageSize[d] = qCeil(10*fieldOfView[d]/spatialResolution[d]);
		points *= imageSize[d];
	}

	Phantom phantom(fieldOfView);

	std::vector<float> image(points);
	float z;
	for(int n=0; n<points; n++)
	{
		float x = (n%imageSize[0]-imageSize[0]/2)*spatialResolution[0]/10;
		float y = ((n/imageSize[0])%imageSize[1] - imageSize[1]/2)*spatialResolution[1]/10;
		if(dimensions>2)
			z = (n/(imageSize[0]*imageSize[1]) - imageSize[2]/2)*spatialResolution[2]/10;

//		if(!(n%333))
//			qWarning("%d %f %f %f", n, x, y, z);

		if(dimensions==3)
			image[n] = phantom.imageDomainSignal(x,y,z);
		else
			image[n] = phantom.imageDomainSignal(x,y);
	}

	writeOctaveFileFloat("test.txt", image, imageSize);
}

void PhantomTest::testFourier_data()
{
	QTest::addColumn<std::vector<float> >("fieldOfView");
	QTest::addColumn<std::vector<float> >("spatialResolution");

	std::vector<float> fieldOfView(2);
	std::vector<float> spatialResolution(2);

	for(int d=0; d<2; d++)
	{
		fieldOfView[d] = 28;
		spatialResolution[d] = 2;
	}

	QTest::newRow("2D") << fieldOfView << spatialResolution;
}

void PhantomTest::testFourier()
{
	QFETCH(std::vector<float>, fieldOfView);
	QFETCH(std::vector<float>, spatialResolution);

	int dimensions = fieldOfView.size();
	std::vector<int> imageSize(dimensions);
	int points = 1;
	for(int d=0; d<dimensions; d++)
	{
		imageSize[d] = qCeil(10*fieldOfView[d]/spatialResolution[d]);
		points *= imageSize[d];
	}

	Phantom phantom(fieldOfView);

	std::vector<complexFloat> kSpace(points);
	float kz;
	for(int n=0; n<points; n++)
	{
		float kx = 5*(n%imageSize[0]-imageSize[0]/2)/(spatialResolution[0]*imageSize[0]);
		float ky = 5*((n/imageSize[0])%imageSize[1] - imageSize[1]/2)/(spatialResolution[1]*imageSize[1]);
		if(dimensions>2)
			kz = (n/(imageSize[0]*imageSize[1]) - imageSize[2]/2)/(spatialResolution[2]*imageSize[2]);

		if(dimensions==3)
			kSpace[n] = phantom.fourierDomainSignal(kx,ky,kz);
		else
			kSpace[n] = phantom.fourierDomainSignal(kx,ky);
	}

	writeOctaveFileFloat("phantomFourierTest", kSpace, imageSize);
}

QTEST_MAIN(PhantomTest)
