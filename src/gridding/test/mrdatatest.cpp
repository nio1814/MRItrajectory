#include "mrdatatest.h"

#include "mrdata.h"

#include <QtTest/QtTest>

void MRdataTest::testFFTshift()
{
	QVector<int> dimensions = QVector<int>() << 4 << 4;

	QVector<complexFloat> signal(16);
	signal[0] = 1;
	signal[1] = 1;
	signal[4] = 1;
	signal[5] = 1;

	QVector<complexFloat> signalShifted(16);
	signalShifted[10] = 1;
	signalShifted[11] = 1;
	signalShifted[14] = 1;
	signalShifted[15] = 1;

	MRdata data(dimensions.toStdVector(), dimensions.size(), signal.toStdVector());

	data.writeToOctave("preshift");
	data.fftShift();
	data.writeToOctave("postshift");
	QVector<complexFloat> dataSignal = QVector<complexFloat>::fromStdVector(data.signal());
	QCOMPARE(dataSignal, signalShifted);
}

QTEST_APPLESS_MAIN(MRdataTest)
