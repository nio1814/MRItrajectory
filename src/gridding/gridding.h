#ifndef GRIDDING_H
#define GRIDDING_H


class Gridding
{
public:
	Gridding();
private:
	std::vector<int> m_imageDimensions;
	float m_oversamplingFactor;
};

#endif // GRIDDING_H
