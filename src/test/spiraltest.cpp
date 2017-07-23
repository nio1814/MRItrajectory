#include "spiraltest.h"

#include <QtTest/QtTest>

extern "C" {
#include "spiral.h"
}

#include "phantom.h"
#include "gridding.h"

void SpiralTest::testGenerate_data()
{
	QTest::addColumn<double>("fieldOfView");
	QTest::addColumn<double>("spatialResolution");
	QTest::addColumn<double>("duration");
	QTest::addColumn<double>("samplingInterval");
	QTest::addColumn<double>("maxGradient");
	QTest::addColumn<double>("maxSlew");
	QTest::addColumn<QVector<QPointF> >("fieldOfViewScale");

	double spatialResolution = 2;
	QVector<QPointF> fieldOfViewScale;

	QTest::newRow("Default") << 28.0 << spatialResolution << 5e-3 << 4e-6 << 4.0 << 15000.0 << fieldOfViewScale;

	fieldOfViewScale.append(QPointF(0,1));

	fieldOfViewScale.append(QPointF(5/spatialResolution, .25));

	QTest::newRow("Variable Density") << 28.0 << spatialResolution << 5e-3 << 4e-6 << 4.0 << 15000.0 << fieldOfViewScale;
}

void SpiralTest::testGenerate()
{
	QFETCH(double, fieldOfView);
	QFETCH(double, spatialResolution);
	QFETCH(double, duration);
	QFETCH(double, samplingInterval);
	QFETCH(double, maxGradient);
	QFETCH(double, maxSlew);
	QFETCH(QVector<QPointF>, fieldOfViewScale);

	VariableDensity* variableDensity = NULL;
	if(fieldOfViewScale.size())
	{
		variableDensity = newVariableDensity();
		for(int n=0; n<fieldOfViewScale.size(); n++)
		{
			QPointF fieldOfViewScalePoint = fieldOfViewScale.at(n);
			addLinearVariableDensityStep(variableDensity, fieldOfViewScalePoint.x(), fieldOfViewScalePoint.y());
		}
	}
	Trajectory* spiral = generateSpirals(variableDensity, fieldOfView, spatialResolution, duration, samplingInterval, 0, Archimedean , 0, fieldOfView, maxGradient, maxSlew);

	saveTrajectory("spiral.trj", spiral);

	qWarning() << "Verifying parameters";
	float threshold = .01;
	for(int d=0; d<2; d++)
	{
		QCOMPARE(spiral->fieldOfView[d], fieldOfView);
		QVERIFY(qAbs(spiral->spatialResolution[d]-spatialResolution) < threshold);
	}
	QCOMPARE(spiral->bases, 1);

	float maxReadoutGradient = qMin(1/(fieldOfView*samplingInterval*4257), maxGradient);
	int waveforms = spiral->storage==StoreAll ? spiral->readouts : spiral->bases;

	QStringList readoutGradientErrorMessages;
	bool readoutGradientExceeded = false;

	char message[128];
	float krMax = 0;
	for(int r=0; r<waveforms; r++)
	{
		for(int n=0; n<spiral->waveformPoints; n++)
		{
			float gradientMagnitude = 0;
			for(int d=0; d<spiral->dimensions; d++)
			{
				float gradientValue = spiral->gradientWaveforms[(d+2*r)*spiral->waveformPoints+n];
				gradientMagnitude += gradientValue*gradientValue;
			}
			gradientMagnitude = qSqrt(gradientMagnitude);
			if(n<spiral->readoutPoints)
			{
				sprintf(message, "|g[%d]| %f limit %f\n", n, gradientMagnitude, maxReadoutGradient);
//				QVERIFY2(gradientMagnitude < maxReadoutGradient+.01, message);
	if(gradientMagnitude>maxReadoutGradient)
				{
						readoutGradientErrorMessages.append(QString("|g[%1,%2]| %3 limit %4\n").arg(r).arg(n).arg(gradientMagnitude).arg(maxReadoutGradient));
					readoutGradientExceeded = true;
				}
				float kr = 0;
				for(int d=0; d<2; d++)
				{
					float k = spiral->kSpaceCoordinates[(d+2*r)+n];
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

	if(readoutGradientExceeded)
	{
		foreach(QString message, readoutGradientErrorMessages)
			qWarning() << message;
	}


	sprintf(message, "max|k| %f expected %f", krMax, 5/spiral->spatialResolution[0]);
	QVERIFY2(krMax>=5/spiral->spatialResolution[0], message);
}

void SpiralTest::testPhantom()
{
	std::vector<float> fieldOfView;
	for(int d=0; d<2; d++)
	{
		fieldOfView.push_back(28);
	}
	Trajectory* spiral = generateSpirals(NULL, fieldOfView[0], 2, 5e-3, 4e-6, 0, Archimedean , 0, fieldOfView[0], 4, 15000);

	std::vector<int> acquisitionSize;
	acquisitionSize.push_back(spiral->readoutPoints);
	acquisitionSize.push_back(spiral->readouts);

	Phantom phantom(fieldOfView);
	MRdata kSpaceData(acquisitionSize, 2);
	for(int n=0; n<spiral->readoutPoints; n++)
	{
		for(int r=0; r<spiral->readouts; r++)
		{
			float k[2];
			trajectoryCoordinates(n, r, spiral, k, NULL);
			int m = spiral->readoutPoints*r + n;
			kSpaceData.setSignalValue(m, phantom.fourierDomainSignal(k[0], k[1]));
//			kSpaceData.setSignalValue(m, 1);
		}
	}

//	saveTrajectory("spiral.trj", spiral);

	Gridding gridding(spiral);
	MRdata* image = gridding.kSpaceToImage(kSpaceData)[0];
//	MRdata* image = gridding.conjugatePhaseForward(kSpaceData);
	image->writeToOctave("spiral.txt");
}

QTEST_MAIN(SpiralTest)
