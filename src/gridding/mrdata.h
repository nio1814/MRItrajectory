#ifndef MRDATA_H
#define MRDATA_H

#include <vector>
#include <complex>

#include <fftw3.h>

typedef std::complex<float> complexFloat;

struct FFTplan
{
	fftwf_plan plan;
	bool inPlace;
};

class MRdata
{
public:
	MRdata(std::vector<int> dimensions, int numImagingDimensions, const std::vector<complexFloat> &data=std::vector<complexFloat>());
	long points() const;

	void setSignalValue(int index, const complexFloat& value);
	complexFloat signalValue(int index) const;
	complexFloat* signal();

	void fftShift();
	FFTplan planFFT(int direction, std::vector<int> outputDimensions);
	void fft(int direction, std::vector<int> outputDimensions=std::vector<int>());
	void crop(std::vector<int> newSize);

	bool writeToOctave(std::string filename);
private:
	std::vector<int> m_dimensions;
	int m_numImagingDimensions;
	std::vector<complexFloat> m_signal;
};

#endif // MRDATA_H
