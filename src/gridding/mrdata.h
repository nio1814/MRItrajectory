#ifndef MRDATA_H
#define MRDATA_H

#include <vector>
#include <complex>

typedef std::complex<float> complexFloat;

class MRdata
{
public:
	MRdata(std::vector<int> dimensions, int numImagingDimensions, const std::vector<complexFloat> &data=std::vector<complexFloat>());
	long points() const;

	void setSignalValue(int index, const complexFloat& value);
	complexFloat signalValue(int index) const;
	complexFloat* signal();

	void fftShift();
	void crop(std::vector<int> newSize);
private:
	int m_numImagingDimensions;
	std::vector<int> m_dimensions;
	std::vector<complexFloat> m_signal;
};

#endif // MRDATA_H
