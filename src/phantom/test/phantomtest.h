#ifndef PHANTOMTEST_H
#define PHANTOMTEST_H

#include <QObject>

class PhantomTest : public QObject
{
	Q_OBJECT
private slots:
	void test3D_data();
	void test3D();
};

#endif // PHANTOMTEST_H
