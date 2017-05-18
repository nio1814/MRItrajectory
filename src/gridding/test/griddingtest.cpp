#include "griddingtest.h"

#include "mrdata.h"
#include "trajectory.h"
#include "gridding.h"

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

	int readoutPoints = randomInteger(10,50);
	qWarning() << "Readout Points:" << readoutPoints;
	int readouts = randomInteger(5,50);
	qWarning() << "Readouts:" << readouts;

	Trajectory trajectory;
	trajectory.dimensions = fieldOfView.size();
	float minResolution = INFINITY;
	for (int d=0; d<trajectory.dimensions; d++)
		minResolution = qMin(minResolution, spatialResolution[d]);

	QVector<float> kSpaceExtent(trajectory.dimensions);
	for (int d=0; d<trajectory.dimensions; d++)
	{
		kSpaceExtent[d] = .5*minResolution/spatialResolution[d];
		trajectory.imageDimensions[d] = 10*fieldOfView[d]/spatialResolution[d];
	}

	QVector<float> kx;
	QVector<float> ky;
	QVector<float> dcf;
	QVector<complexFloat> signal;
	for(int r=0; r<readouts; r++)
	{
		for (int n=0; n<readoutPoints; n++)
		{
			kx.append(randomNumber(-kSpaceExtent[0], kSpaceExtent[0]));
			ky.append(randomNumber(-kSpaceExtent[1], kSpaceExtent[1]));
			dcf.append(randomNumber(0, 1));
			signal.append(complexFloat(randomNumber(-1, 1), randomNumber(-1, 1)));
		}
	}

	QVector<int> trajectorySize = QVector<int>() << trajectory.readoutPoints << trajectory.readouts;
	MRdata kSpaceData(trajectorySize.toStdVector(), trajectory.dimensions, signal.toStdVector());

	Gridding gridding(&trajectory);
	MRdata *griddedKspaceData = gridding.grid(kSpaceData);

}

QTEST_APPLESS_MAIN(GriddingTest)
