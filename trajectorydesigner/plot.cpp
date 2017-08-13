#include "plot.h"

#include <qwt_plot_curve.h>

Plot::Plot(int curves, QWidget *parent) : QwtPlot(parent)
{
	for(int n=0; n<curves; n++)
	{
		m_plotCurve.append(new QwtPlotCurve);
		m_plotCurve.last()->setPen(m_colorOrder[n]);
		m_plotCurve.last()->attach(this);
	}
	setCanvasBackground(Qt::black);
}

void Plot::setYLimits(float yMin, float yMax)
{
	setAxisScale(QwtPlot::yLeft, yMin, yMax);
}

