#ifndef MRDATATEST_H
#define MRDATATEST_H

#include <QObject>

class MRdataTest : public QObject
{
	Q_OBJECT
private slots:
	void testFFTshift();
	void testCrop();
};

#endif // MRDATATEST_H
