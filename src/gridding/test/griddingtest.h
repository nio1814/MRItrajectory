#ifndef GRIDDINGTEST_H
#define GRIDDINGTEST_H

#include <QObject>

class GriddingTest : public QObject
{
	Q_OBJECT
private slots:
	void testForward_data();
	void testForward();
};

#endif // GRIDDINGTEST_H
