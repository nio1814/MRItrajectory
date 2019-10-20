/***************************************************************************

 Copyright (c) 2014 The Board of Trustees of the Leland Stanford Junior University.
 All rights reserved.
 Contact: Okai Addy <noaddy@alumni.stanford.edu>

 This source code is under a BSD 3-Clause License.
 See LICENSE for more information.

To distribute this file, substitute the full license for the above reference.

**************************************************************************/
#include "griddingtest.h"

#include "mrdata.h"
#include "gridding.h"

extern "C"
{
#include "trajectory.h"
}

#include <QtTest/QtTest>

Q_DECLARE_METATYPE(complexFloat)

/*!
 * \brief Generate a random number
 * \param minLimit	Minimum possible number
 * \param maxLimit	Maximum possible number
 */
float randomNumber(float minLimit, float maxLimit)
{
	float range = maxLimit - minLimit;
	return (rand()*range)/RAND_MAX + minLimit;
}

int randomInteger(int minLimit, int maxLimit)
{
	int range = maxLimit - minLimit;
	return rand()%range + minLimit;
}

float meanSquaredError(std::vector<complexFloat> data, std::vector<complexFloat> actual, bool normalize=true)
{
	float energyDifference = 0;
	float energyActual = 0;
	for(size_t n=0; n<actual.size(); n++)
	{
		complexFloat difference = data[n]-actual[n];

		energyDifference += (difference*conj(difference)).real();
		energyActual += (actual[n]*conj(actual[n])).real();
	}

	float mse = std::sqrt(energyDifference);
	if(normalize)
		mse /= std::sqrt(energyActual);

	return mse;
}

void GriddingTest::kernelTest_data()
{
	QTest::addColumn<QVector<float> >("fieldOfView");
	QTest::addColumn<QVector<float> >("spatialResolution");
	QTest::addColumn<float>("oversamplingRatio");

	QVector<float> fieldOfView(2);
	QVector<float> spatialResolution(2);

	fieldOfView = QVector<float>() << randomNumber(10, 48);
	fieldOfView << fieldOfView;
  spatialResolution = QVector<float>() << randomNumber(.8f, 4.0f);
	spatialResolution << spatialResolution;

	QTest::newRow("Isotropic") << fieldOfView << spatialResolution << randomNumber(1.0, 3.0);
}

void GriddingTest::kernelTest()
{
	QFETCH(QVector<float>, fieldOfView);
	QFETCH(QVector<float>, spatialResolution);
	QFETCH(float, oversamplingRatio);

	Trajectory trajectory;
  trajectory.numDimensions = fieldOfView.size();
  for(int d=0; d<trajectory.numDimensions; d++)
	{
		trajectory.spatialResolution[d] = spatialResolution[d];
		trajectory.fieldOfView[d] = fieldOfView[d];
	}

	Gridding gridding(&trajectory, oversamplingRatio);

	qInfo() << "Lookup table";
	qInfo() << QVector<float>::fromStdVector(gridding.m_kernelLookupTable).toList();
}

void GriddingTest::testForward_data()
{
	QTest::addColumn<QVector<float> >("fieldOfView");
	QTest::addColumn<QVector<float> >("spatialResolution");
	QTest::addColumn<int>("readouts");
	QTest::addColumn<int>("readoutPoints");
	QTest::addColumn<QVector<float> >("kSpaceCoordinates");
	QTest::addColumn<QVector<float> >("densityCompensation");
	QTest::addColumn<QVector<complexFloat> >("signal");
	QTest::addColumn<float>("oversamplingRatio");

	QVector<float> fieldOfView(2);
	QVector<float> spatialResolution(2);

	fieldOfView = QVector<float>() << randomNumber(10, 48);
	fieldOfView << fieldOfView;
  spatialResolution = QVector<float>() << randomNumber(.8f, 4.0f);
	spatialResolution << spatialResolution;

	float kSpaceExtent[3];
  int numDimensions = 2;
  for (int d=0; d<numDimensions; d++)
		kSpaceExtent[d] = 5/spatialResolution[d];

	QVector<float> kSpaceCoordinates;
	QVector<float> densityCompensation;
	QVector<complexFloat> signal;

	int readoutPoints = randomInteger(10,50);
	int readouts = randomInteger(5,50);
	for(int r=0; r<readouts; r++)
	{
		for (int n=0; n<readoutPoints; n++)
		{
      for(int d=0; d<numDimensions; d++)
			{
				kSpaceCoordinates.append(randomNumber(-kSpaceExtent[d], kSpaceExtent[d]));
			}
			densityCompensation.append(randomNumber(0,1));
			signal.append(complexFloat(randomNumber(-1, 1), randomNumber(-1, 1)));
		}
	}

	QTest::newRow("Isotropic Random") << fieldOfView << spatialResolution << readoutPoints << readouts << kSpaceCoordinates << densityCompensation << signal << randomNumber(1.0, 3.0);

	for(int d=0; d<2; d++)
	{
		fieldOfView[d] = 28;
		spatialResolution[d] = 2;
	}

	signal.clear();
	signal.append(complexFloat(1,0));

	kSpaceCoordinates.clear();
	kSpaceCoordinates.append(0);
	kSpaceCoordinates.append(0);

	densityCompensation.clear();
	densityCompensation.append(1);

	QTest::newRow("DC") << fieldOfView << spatialResolution << 1 << 1 << kSpaceCoordinates << densityCompensation << signal << 2.0f;

	signal.append(complexFloat(1,0));
	kSpaceCoordinates.append(1);
	kSpaceCoordinates.append(0);
	densityCompensation.append(1);
}

void GriddingTest::testForward()
{
	QFETCH(QVector<float>, fieldOfView);
	QFETCH(QVector<float>, spatialResolution);
	QFETCH(int, readouts);
	QFETCH(int, readoutPoints);
	QFETCH(QVector<float>, kSpaceCoordinates);
	QFETCH(QVector<float>, densityCompensation);
	QFETCH(QVector<complexFloat>, signal);
	QFETCH(float, oversamplingRatio);

	qWarning() << "Oversampling Ratio:" << oversamplingRatio;

	Trajectory trajectory;
  trajectory.numReadoutPoints = readoutPoints;
  qWarning() << "Readout Points:" << trajectory.numReadoutPoints;
  trajectory.numReadouts = readouts;
  qWarning() << "Readouts:" << trajectory.numReadouts;
  trajectory.numDimensions = fieldOfView.size();
  allocateTrajectory(&trajectory, trajectory.numReadoutPoints, 0, trajectory.numDimensions, trajectory.numReadouts, trajectory.numReadouts, STORE_ALL);
	float minResolution = INFINITY;
  for (int d=0; d<trajectory.numDimensions; d++)
	{
		trajectory.spatialResolution[d] = spatialResolution[d];
		trajectory.fieldOfView[d] = fieldOfView[d];
		minResolution = qMin(minResolution, spatialResolution[d]);
		trajectory.imageDimensions[d] = 10*fieldOfView[d]/spatialResolution[d];
	}

  for(int r=0; r<trajectory.numReadouts; r++)
	{
    for (int n=0; n<trajectory.numReadoutPoints; n++)
		{
			int m = r*readoutPoints+n;
			float k[2] = {kSpaceCoordinates[2*m], kSpaceCoordinates[2*m+1]};
			setTrajectoryPoint(n, r, &trajectory, k, randomNumber(0, 1));
		}
	}

  QVector<int> trajectorySize = QVector<int>() << trajectory.numReadoutPoints << trajectory.numReadouts;
  MRdata kSpaceData(trajectorySize.toStdVector(), trajectory.numDimensions, signal.toStdVector());

	kSpaceData.writeToOctave("kspace.txt");

	Gridding gridding(&trajectory, oversamplingRatio);
	std::vector<MRdata*> griddedData = gridding.kSpaceToImage(kSpaceData, true);

	MRdata* griddedKspace = griddedData[0];
	griddedKspace->writeToOctave("griddedkspace.txt");

	MRdata* image = griddedData[1];
	image->writeToOctave("grid.txt");

	MRdata* conjugatePhaseImage = gridding.conjugatePhaseForward(kSpaceData);

	conjugatePhaseImage->writeToOctave("conjugatephase.txt");

	float error = meanSquaredError(image->signal(), conjugatePhaseImage->signal());

	QString message = QString("error %1").arg(error);

	qInfo() << "Image error " << error;
	QVERIFY(error<.04);

	MRdata* kSpaceDataForwardInverse = gridding.imageToKspace(*image);

	QCOMPARE(kSpaceData.points(), kSpaceDataForwardInverse->points());

	error = meanSquaredError(kSpaceDataForwardInverse->signal(), kSpaceData.signal());

	message = QString("error %1").arg(error);
	QVERIFY2(error<.04, message.toStdString().c_str());
}

QTEST_APPLESS_MAIN(GriddingTest)
