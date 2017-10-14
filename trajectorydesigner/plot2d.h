#ifndef PLOT2D_H
#define PLOT2D_H

#include "plot.h"

class Plot2D : public Plot
{
public:
	Plot2D();

	void setSamples(QVector<QPointF> samples);
};

#endif // PLOT2D_H
