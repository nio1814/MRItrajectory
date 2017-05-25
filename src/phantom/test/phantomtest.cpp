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
	stream << "# ndims: 3\n";
	QString line = QString(" %1 %2 %3\n").arg(imageSize[0]).arg(imageSize[1]).arg(imageSize[2]);
	stream << line;

	for (size_t n=0; n<data.size(); n++)
		stream << " " << data[n] << endl;

	file.close();
}

void PhantomTest::test3D_data()
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
	QTest::newRow("Default") << fieldOfView << spatialResolution;
}

void PhantomTest::test3D()
{
	QFETCH(std::vector<float>, fieldOfView);
	QFETCH(std::vector<float>, spatialResolution);

	std::vector<int> imageSize(3);
	int points = 1;
	for(int d=0; d<3; d++)
	{
		imageSize[d] = qCeil(10*fieldOfView[d]/spatialResolution[d]);
		points *= imageSize[d];
	}

	Phantom phantom(fieldOfView);

	std::vector<float> image(points);
	for(int n=0; n<points; n++)
	{
		float x = (n%imageSize[0]-imageSize[0]/2)*spatialResolution[0]/10;
		float y = ((n/imageSize[0])%imageSize[1] - imageSize[1]/2)*spatialResolution[1]/10;
		float z = (n/(imageSize[0]*imageSize[1]) - imageSize[2]/2)*spatialResolution[2]/10;

		if(!(n%333))
			qWarning("%d %f %f %f", n, x, y, z);

		image[n] = phantom.imageDomainSignal(x,y,z);
	}

	writeOctaveFileFloat("test.txt", image, imageSize);
}

QTEST_MAIN(PhantomTest)
