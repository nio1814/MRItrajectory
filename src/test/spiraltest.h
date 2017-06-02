#ifndef SPIRALTEST_H
#define SPIRALTEST_H

#include <QObject>

class SpiralTest : public QObject
{
	Q_OBJECT
private slots:
	void testPhantom();
	void testGenerate_data();
	void testGenerate();
};

#endif // SPIRALTEST_H
