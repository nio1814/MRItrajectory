#ifndef CONESTEST_H
#define CONESTEST_H

#include <QObject>

class ConesTest : public QObject
{
	Q_OBJECT
private slots:
	void generateTest_data();
	void generateTest();
};

#endif // CONESTEST_H
