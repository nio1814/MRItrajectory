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
		energyDifference += abs(data[n]-actual[n]);
		energyActual += (actual[n]*conj(actual[n])).real();
	}

	float mse = std::sqrt(energyDifference);
	if(normalize)
		mse /= std::sqrt(energyActual);

	return mse;
}

void GriddingTest::testForward_data()
{
	QTest::addColumn<QVector<float> >("fieldOfView");
	QTest::addColumn<QVector<float> >("spatialResolution");
	QTest::addColumn<float>("oversamplingRatio");

	QVector<float> fieldOfView(2);
	QVector<float> spatialResolution(2);

	fieldOfView = QVector<float>() << randomNumber(10, 48);
	fieldOfView << fieldOfView;
	spatialResolution = QVector<float>() << randomNumber(.8, 4.0);
	spatialResolution << spatialResolution;

	QTest::newRow("Isotropic") << fieldOfView << spatialResolution << randomNumber(1.0, 3.0);
}

void GriddingTest::testForward()
{
	QFETCH(QVector<float>, fieldOfView);
	QFETCH(QVector<float>, spatialResolution);
	QFETCH(float, oversamplingRatio);

	qWarning() << "Oversampling Ratio:" << oversamplingRatio;

	Trajectory trajectory;
	trajectory.readoutPoints = randomInteger(10,50);
	qWarning() << "Readout Points:" << trajectory.readoutPoints;
	trajectory.readouts = randomInteger(5,50);
	qWarning() << "Readouts:" << trajectory.readouts;
	trajectory.dimensions = fieldOfView.size();
	allocateTrajectory(&trajectory, trajectory.readoutPoints, 0, trajectory.dimensions, trajectory.readouts, trajectory.readouts, StoreAll);
	float minResolution = INFINITY;
	for (int d=0; d<trajectory.dimensions; d++)
	{
		trajectory.spatialResolution[d] = spatialResolution[d];
		trajectory.fieldOfView[d] = fieldOfView[d];
		minResolution = qMin(minResolution, spatialResolution[d]);
	}

	QVector<float> kSpaceExtent(trajectory.dimensions);
	for (int d=0; d<trajectory.dimensions; d++)
	{
		kSpaceExtent[d] = 5/spatialResolution[d];
		trajectory.imageDimensions[d] = 10*fieldOfView[d]/spatialResolution[d];
	}

//	QVector<float> kx;
//	QVector<float> ky;
	QVector<float> k(trajectory.dimensions);
//	QVector<float> dcf;
	QVector<complexFloat> signal;
	for(int r=0; r<trajectory.readouts; r++)
	{
		for (int n=0; n<trajectory.readoutPoints; n++)
		{
			k.clear();
			for(int d=0; d<trajectory.dimensions; d++)
			{
				k.append(randomNumber(-kSpaceExtent[d], kSpaceExtent[d]));
			}
			setTrajectoryPoint(n, r, &trajectory, k.data(), randomNumber(0, 1));
			signal.append(complexFloat(randomNumber(-1, 1), randomNumber(-1, 1)));
		}
	}

	QVector<int> trajectorySize = QVector<int>() << trajectory.readoutPoints << trajectory.readouts;
	MRdata kSpaceData(trajectorySize.toStdVector(), trajectory.dimensions, signal.toStdVector());

	Gridding gridding(&trajectory);
	MRdata* image = gridding.kSpaceToImage(kSpaceData);

	image->writeToOctave("grid.txt");

	MRdata* conjugatePhaseImage = gridding.conjugatePhaseForward(kSpaceData);

	conjugatePhaseImage->writeToOctave("conjugatephase.txt");

	float error = meanSquaredError(image->signal(), conjugatePhaseImage->signal());

	QString message = QString("error %1").arg(error);

	QVERIFY2(error<.04, message.toStdString().c_str());

	MRdata* kSpaceDataForwardInverse = gridding.imageToKspace(*image);

	QCOMPARE(kSpaceData.points(), kSpaceDataForwardInverse->points());

	error = meanSquaredError(kSpaceDataForwardInverse->signal(), kSpaceData.signal());

	message = QString("error %1").arg(error);
	QVERIFY2(error<.04, message.toStdString().c_str());
}

QTEST_APPLESS_MAIN(GriddingTest)
