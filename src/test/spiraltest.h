#ifndef SPIRALTEST_H
#define SPIRALTEST_H

#include <QObject>

class SpiralTest : public QObject
{
	Q_OBJECT
private slots:
	void testGenerate_data();
	void testGenerate();
	void testPhantom();
};

#endif // SPIRALTEST_H
