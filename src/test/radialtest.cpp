#include "radialtest.h"

extern "C" {
#include "radial.h"
#include "trajectory.h"
}

#include "mrdata.h"
#include "gridding.h"
#include "phantom.h"

#include <QtTest/QtTest>

void RadialTest::testGenerate_data()
{
	QTest::addColumn<QVector<float> >("fieldOfView");
	QTest::addColumn<QVector<float> >("spatialResolution");
	QTest::addColumn<double>("duration");
	QTest::addColumn<double>("samplingInterval");
	QTest::addColumn<double>("maxGradient");
	QTest::addColumn<double>("maxSlew");

	QVector<float> fieldOfView = QVector<float>() << 28 << 28;
	QVector<float> spatialResolution = QVector<float>() << 2 << 2;

	QTest::newRow("Default") << fieldOfView << spatialResolution << 5e-3 << 4e-6 << 4.0 << 15000.0;
}

void RadialTest::testGenerate()
{
	QFETCH(QVector<float>, fieldOfView);
	QFETCH(QVector<float>, spatialResolution);
	QFETCH(double, duration);
	QFETCH(double, samplingInterval);
	QFETCH(double, maxGradient);
	QFETCH(double, maxSlew);

	Trajectory* radial = generateRadial(fieldOfView[0], fieldOfView[1], InverseEllipticalShape, spatialResolution[0], spatialResolution[1], InverseEllipticalShape, 1, 1, maxGradient, maxSlew, samplingInterval);

	saveTrajectory("radial.trj", radial);

	qWarning() << "Verifying parameters";
	float threshold = .01;
	for(int d=0; d<2; d++)
	{
		QCOMPARE(radial->fieldOfView[d], fieldOfView[d]);
		QVERIFY(qAbs(radial->spatialResolution[d]-spatialResolution[d]) < threshold);
	}
	QCOMPARE(radial->bases, 1);

	float fieldOfViewMax = fmax(fieldOfView[0], fieldOfView[1]);
	float maxReadoutGradient = qMin(1/(fieldOfViewMax*samplingInterval*4257), maxGradient);
	int waveforms;

	char message[128];
	float krMax = 0;
	for(int r=0; r<waveforms; r++)
	{
		for(int n=0; n<radial->waveformPoints; n++)
		{
			float gradientMagnitude = 0;
			for(int d=0; d<radial->dimensions; d++)
			{
				float gradientValue = radial->gradientWaveforms[(d+2*r)*radial->waveformPoints+n];
				gradientMagnitude += gradientValue*gradientValue;
			}
			gradientMagnitude = qSqrt(gradientMagnitude);
			if(n<radial->readoutPoints)
			{
				sprintf(message, "|g[%d]| %f limit %f\n", n, gradientMagnitude, maxReadoutGradient);
				QVERIFY2(gradientMagnitude < maxReadoutGradient+.01, message);
				float kr = 0;
				for(int d=0; d<2; d++)
				{
					float k = radial->kSpaceCoordinates[(d+2*r)+n];
					kr += k*k;
				}
				kr = qSqrt(kr);
				krMax = qMax(kr, krMax);
			}
			else
			{
				sprintf(message, "|g[%d,%d]| %f limit %f\n", r, n, gradientMagnitude, maxGradient);
				QVERIFY2(gradientMagnitude < maxGradient, message);
			}
		}
	}

	sprintf(message, "max|k| %f expected %f", krMax, 5/radial->spatialResolution[0]);
	QVERIFY2(krMax>=5/radial->spatialResolution[0], message);
}

void RadialTest::testPhantom()
{
	std::vector<float> fieldOfView;
	std::vector<float> spatialResolution;
	for(int d=0; d<2; d++)
	{
		fieldOfView.push_back(28);
		spatialResolution.push_back(2);
	}

	Trajectory* radial = generateRadial(fieldOfView[0], fieldOfView[1], InverseEllipticalShape, spatialResolution[0], spatialResolution[1], InverseEllipticalShape, 1, 1, 4, 15000, 4e-6);

	std::vector<int> acquisitionSize;
	acquisitionSize.push_back(radial->readoutPoints);
	acquisitionSize.push_back(radial->readouts);

	Phantom phantom(fieldOfView);
	MRdata kSpaceData(acquisitionSize, 2);
	for(int n=0; n<radial->readoutPoints; n++)
	{
		for(int r=0; r<radial->readouts; r++)
		{
			float k[2];
			trajectoryCoordinates(n, r, radial, k, NULL);
			int m = radial->readoutPoints*r + n;
//			if(n<69)
			kSpaceData.setSignalValue(m, phantom.fourierDomainSignal(k[0], k[1]));
//			kSpaceData.setSignalValue(m, 1);
		}
	}

	saveTrajectory("radial.trj", radial);

	Gridding gridding(radial);
	MRdata* image = gridding.kSpaceToImage(kSpaceData)[0];
//	MRdata* image = gridding.conjugatePhaseForward(kSpaceData);
	image->writeToOctave("radial.txt");
}


QTEST_MAIN(RadialTest)
