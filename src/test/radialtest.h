#ifndef RADIALTEST_H
#define RADIALTEST_H

#include <QObject>

class RadialTest : public QObject
{
	Q_OBJECT
private slots:
	void testGenerate_data();
	void testGenerate();
	void testPhantom();
};

#endif // RADIALTEST_H
