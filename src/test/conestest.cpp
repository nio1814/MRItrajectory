#include "conestest.h"

extern "C" {
#include "cones.h"
}

#include <QtTest/QtTest>

Q_DECLARE_METATYPE(InterConeCompensation)
Q_DECLARE_METATYPE(WaveformStorageType)

void ConesTest::generateTest_data()
{
	QTest::addColumn<QVector<float> >("fieldOfView");
	QTest::addColumn<QVector<float> >("spatialResolution");
	QTest::addColumn<int>("bases");
	QTest::addColumn<int>("rotatable");
		QTest::addColumn<InterConeCompensation>("interconeCompensationType");
	QTest::addColumn<double>("duration");
	QTest::addColumn<double>("samplingInterval");
	QTest::addColumn<double>("maxGradient");
	QTest::addColumn<double>("maxSlew");
	QTest::addColumn<WaveformStorageType>("storage");

	QTest::newRow("Isotropic") << (QVector<float>() << 28 << 28) << (QVector<float>() << 2 << 2) << 16 << 1 << NoCompensation << 8.4e-3 << 4e-6 << 4.0 << 15000.0 << StoreBasis;

//	QTest::newRow("Anisotropic field of view") << (QVector<float>() << 28 << 14 ) << (QVector<float>() << 2 << 2) << 16 << 1 << NoCompensation << 5e-3 << 4e-6 << 4.0 << 15000.0 << StoreAll;

//	QTest::newRow("Cones") << (QVector<float>() << 28 << 14) << (QVector<float>() << 1.2 << 1.25) << 32 << 1 << 1 << 2.8e-3 << 4e-6 << 4.0 << 15000.0 << StoreAll;
}

void ConesTest::generateTest()
{
	QFETCH(QVector<float>, fieldOfView);
	QFETCH(QVector<float>, spatialResolution);
	QFETCH(int, bases);
	QFETCH(int, rotatable);
	QFETCH(InterConeCompensation, interconeCompensationType);
	QFETCH(double, duration);
	QFETCH(double, samplingInterval);
	QFETCH(double, maxGradient);
	QFETCH(double, maxSlew);
	QFETCH(WaveformStorageType, storage);

	float filterFieldOfView = qMax(fieldOfView[0], fieldOfView[1]);

	struct Cones* cones = generateCones(fieldOfView[0], fieldOfView[1], NULL, spatialResolution[0], spatialResolution[1], bases, rotatable, interconeCompensationType, duration, samplingInterval, filterFieldOfView, maxGradient, maxSlew, storage);
	struct Trajectory *trajectory = &cones->trajectory;
	saveTrajectory("cones.trj", trajectory);

	qWarning() << "Verifying parameters";
	QCOMPARE(trajectory->fieldOfView[0], fieldOfView[0]);
	QCOMPARE(trajectory->fieldOfView[1], fieldOfView[0]);
	QCOMPARE(trajectory->fieldOfView[2], fieldOfView[1]);
	float threshold = .01;
	QVERIFY(qAbs(trajectory->spatialResolution[0]-spatialResolution[0]) < threshold);
	QVERIFY(qAbs(trajectory->spatialResolution[1]- spatialResolution[0]) < threshold);
	QVERIFY(qAbs(trajectory->spatialResolution[2]- spatialResolution[1]) < threshold);
	QCOMPARE(trajectory->bases, bases);

	float maxReadoutGradient = 1/(filterFieldOfView*samplingInterval*4257);
	int waveforms;
	if(storage==StoreBasis)
		waveforms = trajectory->bases;
	else
		waveforms = trajectory->readouts;

	char message[128];
	float kxyMax;
	float kzMax;
	struct ConesInterpolation* schedule = &cones->interpolation;
	for(int r=0; r<waveforms; r++)
	{
		int b;
		if(storage=StoreBasis)
			b = r;
		else
			b = schedule->basis[r];
		for(int n=0; n<trajectory->waveformPoints; n++)
		{
			float gradientMagnitude = 0;
			for(int d=0; d<trajectory->dimensions; d++)
			{
				float* gradientWaveform = trajectoryGradientWaveform(trajectory, r, d);
				float gradientValue = gradientWaveform[n];
				gradientMagnitude += gradientValue*gradientValue;
			}
			gradientMagnitude = qSqrt(gradientMagnitude);
			if(n<cones->basisReadoutPoints[b])
			{
				sprintf(message, "|g| %f limit %f\n", gradientMagnitude, maxReadoutGradient);
				QVERIFY2(gradientMagnitude < maxReadoutGradient+.01, message);
				float k[3];
				trajectoryCoordinates(n, b, trajectory, k, NULL);
				float kxy = qSqrt(k[0]*k[0]+k[1]*k[1]);
				kxyMax = qMax(kxy, kxyMax);
				kzMax = qMax(k[2], kzMax);
			}
			else
			{
				sprintf(message, "|g[%d,%d]| %f limit %f\n", r, n, gradientMagnitude, maxGradient);
				QVERIFY2(gradientMagnitude < maxGradient, message);
			}
		}
	}

	sprintf(message, "max|kxy| %f expected %f", kxyMax, 5/trajectory->spatialResolution[0]);
	QVERIFY2(kxyMax>=5/trajectory->spatialResolution[0], message);
	sprintf(message, "max|kz| %f expected %f", kzMax, 5/trajectory->spatialResolution[1]);
	QVERIFY2(kzMax>=5/trajectory->spatialResolution[1], message);
}

QTEST_MAIN(ConesTest)
