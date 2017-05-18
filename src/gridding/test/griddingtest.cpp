#include "griddingtest.h"

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

}

QTEST_APPLESS_MAIN(GriddingTest)