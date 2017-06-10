#ifndef RADIALTEST_H
#define RADIALTEST_H

#include <QObject>

class RadialTest : public QObject
{
	Q_OBJECT
private slots:
	void testGenerate2D_data();
	void testGenerate2D();
	void testPhantom();
};

#endif // RADIALTEST_H
