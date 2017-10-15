#include "variabledensityplot.h"

#include <QPainter>

VariableDensityPlot::VariableDensityPlot(QWidget *parent) : QWidget(parent)
{
	setAutoFillBackground(true);

	QPalette palette;
	palette.setColor(QPalette::Background, Qt::black);
	setPalette(palette);
}

void VariableDensityPlot::paintEvent(QPaintEvent *event)
{
	Q_UNUSED(event);

	QPainter painter(this);

	QPen pen;
	pen.setColor(Qt::white);
	pen.setColor(Qt::yellow);
	pen.setWidth(5);
	painter.setPen(pen);
	QLinearGradient gradient;
	painter.setBrush(gradient);
	QPointF pointPrevious;
	for(int n=0; n<m_points.size(); n++)
	{
		QPainterPath path;
		QPointF pointCurrent = mapToWindow(m_points[n]);
		path.addEllipse(pointCurrent, 5, 5);
		painter.drawPath(path);
		if(n)
			painter.drawLine(pointPrevious, pointCurrent);
		pointPrevious = pointCurrent;
	}
}

QPointF VariableDensityPlot::mapToWindow(QPointF point)
{
	float scale = .9;
	float offsetFactor = (1-scale)*.5;
	float x = (point.x()*scale + offsetFactor)*width();
	float y = (1-point.y() + offsetFactor)*height();

	return QPointF(x,y);
}


