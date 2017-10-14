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
	pen.setColor(Qt::yellow);
	pen.setWidth(10);
	painter.setPen(pen);
	QLinearGradient gradient;
	painter.setBrush(gradient);
	QPointF pointPrevious;
	for(int n=0; n<m_points.size(); n++)
	{
		QPainterPath path;
		float x = m_points[n].x()*width();
		float y = (1-m_points[n].y())*height();
		path.addEllipse(x, y, 5, 5);
		painter.drawPath(path);
		pointPrevious.setX(x);
		pointPrevious.setY(y);
		QPointF pointCurrent(x,y);
		if(n)
			painter.drawLine(pointPrevious, pointCurrent);
	}
}


