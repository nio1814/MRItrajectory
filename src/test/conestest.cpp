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
	QTest::addColumn<QVector<QPointF> >("fieldOfViewScale");
  QTest::addColumn<int>("numReadouts");

	QVector<QPointF> fieldOfViewScale;

  QTest::newRow("Isotropic - basis") << (QVector<float>() << 28 << 28) << (QVector<float>() << 2 << 2) << 16 << 1 << NO_COMPENSATION << 8.4e-3 << 4e-6 << 4.0 << 15000.0 << STORE_BASIS << fieldOfViewScale << 16;

//  QTest::newRow("Isotropic - full") << (QVector<float>() << 28 << 28) << (QVector<float>() << 2 << 2) << 16 << 1 << NoCompensation << 8.4e-3 << 4e-6 << 4.0 << 15000.0 << STORE_ALL << fieldOfViewScale;

//  QTest::newRow("Anisotropic field of view") << (QVector<float>() << 28 << 14 ) << (QVector<float>() << 2 << 2) << 16 << 1 << NoCompensation << 8.4e-3 << 4e-6 << 4.0 << 15000.0 << STORE_ALL << fieldOfViewScale;

//  QTest::newRow("Anisotropic spatial resolution") << (QVector<float>() << 28 << 28 ) << (QVector<float>() << 2 << 4) << 48 << 1 << NoCompensation << 8.4e-3 << 4e-6 << 4.0 << 15000.0 << STORE_ALL << fieldOfViewScale;

//	fieldOfViewScale.append(QPointF(0,1));
//	fieldOfViewScale.append(QPointF(2.5,0.5));

//  QTest::newRow("Variable Density") << (QVector<float>() << 28 << 28 ) << (QVector<float>() << 2 << 2) << 32 << 1 << NoCompensation << 8.4e-3 << 4e-6 << 4.0 << 15000.0 << STORE_ALL << fieldOfViewScale;

  QTest::newRow("Whole Heart Coronary") << (QVector<float>() << 28 << 14) << (QVector<float>() << 1.2 << 1.25) << 16 << 1 << Compensation1 << 2.18e-3 << 4e-6 << 4.0 << 15000.0 << STORE_ALL << fieldOfViewScale << 10405;
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
	QFETCH(QVector<QPointF>, fieldOfViewScale);
  QFETCH(int, numReadouts);

	float filterFieldOfView = qMax(fieldOfView[0], fieldOfView[1]);

  VariableDensity* variableDensity = nullptr;
	if(fieldOfViewScale.size())
	{
		variableDensity = newVariableDensity();
		for(int n=0; n<fieldOfViewScale.size(); n++)
		{
			QPointF fieldOfViewScalePoint = fieldOfViewScale.at(n);
			addLinearVariableDensityStep(variableDensity, fieldOfViewScalePoint.x(), fieldOfViewScalePoint.y());
		}
	}

	struct Cones* cones = generateCones(fieldOfView[0], fieldOfView[1], variableDensity, spatialResolution[0], spatialResolution[1], bases, rotatable, interconeCompensationType, duration, samplingInterval, filterFieldOfView, maxGradient, maxSlew, storage);
  struct Trajectory *trajectory = cones->trajectory;
	saveTrajectory("cones.trj", trajectory);

	qWarning() << "Verifying parameters";
	QCOMPARE(trajectory->fieldOfView[0], fieldOfView[0]);
	QCOMPARE(trajectory->fieldOfView[1], fieldOfView[0]);
	QCOMPARE(trajectory->fieldOfView[2], fieldOfView[1]);
	float threshold = .01;
	QVERIFY(qAbs(trajectory->spatialResolution[0]-spatialResolution[0]) < threshold);
	QVERIFY(qAbs(trajectory->spatialResolution[1]- spatialResolution[0]) < threshold);
	QVERIFY(qAbs(trajectory->spatialResolution[2]- spatialResolution[1]) < threshold);
  QCOMPARE(trajectory->numBases, bases);

	float maxReadoutGradient = 1/(filterFieldOfView*samplingInterval*4257);
  int numWaveformsTotal;
	if(storage==STORE_BASIS)
    numWaveformsTotal = trajectory->numBases;
	else
  {
    numWaveformsTotal = trajectory->numReadouts;
    QVERIFY(abs(numWaveformsTotal - numReadouts) < 10);
  }

	char message[128];
	float kxyMax;
	float kzMax;
  struct ConesInterpolation* schedule = cones->interpolation;
  for(int r=0; r<numWaveformsTotal; r++)
	{
		int b;
		if(storage==STORE_BASIS)
			b = r;
		else
			b = schedule->basis[r];
		for(int n=0; n<trajectory->numWaveformPoints; n++)
		{
			float gradientMagnitude = 0;
			for(int d=0; d<trajectory->numDimensions; d++)
			{
				float* gradientWaveform = trajectoryGradientWaveform(trajectory, r, d);
				float gradientValue = gradientWaveform[n];
				gradientMagnitude += gradientValue*gradientValue;
			}
			gradientMagnitude = qSqrt(gradientMagnitude);
      if(n<cones->numBasisReadoutPoints[b])
			{
				sprintf(message, "|g| %f limit %f\n", gradientMagnitude, maxReadoutGradient);
				QVERIFY2(gradientMagnitude < maxReadoutGradient+.01, message);
				float k[3];
				float densityCompensation;
				trajectoryCoordinates(n, r, trajectory, k, &densityCompensation);
				float kxy = hypotf(k[0],+k[1]);
				kxyMax = qMax(kxy, kxyMax);
				kzMax = qMax(qAbs(k[2]), kzMax);

				sprintf(message, "dcf[%d,%d]=%f\n", r, n, densityCompensation);
				QVERIFY2(densityCompensation>=0, message);
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
	QVERIFY2(kzMax>=5/trajectory->spatialResolution[2], message);
}

QTEST_MAIN(ConesTest)
