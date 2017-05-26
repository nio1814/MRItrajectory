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

void MRdataTest::testCrop()
{
	QVector<int> dimensions = QVector<int>() << 4 << 4;
	QVector<complexFloat> signal;
	for(int n=0; n<16; n++)
		signal.append(n);

	QVector<complexFloat> signalCropped(4);
	signalCropped[0] = signal[5];
	signalCropped[1] = signal[6];
	signalCropped[2] = signal[9];
	signalCropped[3] = signal[10];

	MRdata data(dimensions.toStdVector(), dimensions.size(), signal.toStdVector());

	QVector<int> dimensionsCrop = QVector<int>() << 2 << 2;
	data.writeToOctave("precrop");
	data.crop(dimensionsCrop.toStdVector());
	data.writeToOctave("postcrop");

	QVector<complexFloat> dataSignal = QVector<complexFloat>::fromStdVector(data.signal());
	QCOMPARE(dataSignal, signalCropped);
}

QTEST_APPLESS_MAIN(MRdataTest)
