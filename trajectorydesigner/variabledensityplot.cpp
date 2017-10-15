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

	QPen pen;
	int gridLines = 5;
	for(int n=0; n<gridLines+1; n++)
	{
		if(n)
		{
			pen.setColor(Qt::gray);
			pen.setStyle(Qt::DashLine);
		}
		else
		{
			pen.setColor(Qt::white);
			pen.setStyle(Qt::SolidLine);
		}
		float r = n/(float)gridLines;
		drawLine(QPointF(r,0), QPointF(r,1), pen);
		drawLine(QPointF(0,r), QPointF(1,r), pen);
	}

	QPainter painter(this);

	pen.setColor(Qt::yellow);
	pen.setWidth(5);
	pen.setStyle(Qt::SolidLine);
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
	float y = ((1-point.y())*scale + offsetFactor)*height();

	return QPointF(x,y);
}

void VariableDensityPlot::drawLine(QPointF point1, QPointF point2, const QPen& pen)
{
	point1 = mapToWindow(point1);
	point2 = mapToWindow(point2);

	QPainter painter(this);
	painter.setPen(pen);
	painter.drawLine(point1, point2);
}


