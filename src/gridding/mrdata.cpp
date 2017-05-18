#include "mrdata.h"

MRdata::MRdata(std::vector<int> dimensions, int numImagingDimensions, const std::vector<complexFloat> &data) : m_numImagingDimensions(numImagingDimensions)
{
	m_signal = data;
}

