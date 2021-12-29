#include "plot2d.h"

#include <qwt_plot_grid.h>
#include <qwt_plot_curve.h>

Plot2D::Plot2D(QWidget *parent) : Plot(1, parent)
{
	setCanvasBackground(Qt::black);

	QScopedPointer<QwtPlotGrid> trajectoryPlotGrid(new QwtPlotGrid);
	trajectoryPlotGrid->attach(this);
	trajectoryPlotGrid->setMajorPen(Qt::gray, 0, Qt::DotLine);

	m_plotCurve.append(QSharedPointer<QwtPlotCurve>::create());
	m_plotCurve.last()->setPen(Qt::green);
	m_plotCurve.last()->attach(this);
}

void Plot2D::setSamples(QVector<QPointF> samples)
{
	m_plotCurve[0]->setSamples(samples);
}
