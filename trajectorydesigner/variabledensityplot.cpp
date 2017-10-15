#include "variabledensityplot.h"

#include <QPainter>
#include <QMouseEvent>
#include <QDebug>

extern "C"
{
#include "variabledensity.h"
}

VariableDensityPlot::VariableDensityPlot(QWidget *parent) : QWidget(parent),
	m_variableDensity(newVariableDensity())
{
	setFocus();
	setAutoFillBackground(true);

	QPalette palette;
	palette.setColor(QPalette::Background, Qt::black);
	setPalette(palette);

	addPoint(0, 1);
	addPoint(1, 1);
}

void VariableDensityPlot::paintEvent(QPaintEvent *event)
{
	Q_UNUSED(event);

	QPen penGridLine;
	int gridLines = 5;
	for(int n=0; n<gridLines+1; n++)
	{
		if(n)
		{
			penGridLine.setColor(Qt::gray);
			penGridLine.setStyle(Qt::DashLine);
		}
		else
		{
			penGridLine.setColor(Qt::white);
			penGridLine.setStyle(Qt::SolidLine);
		}
		float r = n/(float)gridLines;
		drawLine(QPointF(r,0), QPointF(r,1), penGridLine);
		drawLine(QPointF(0,r), QPointF(1,r), penGridLine);
	}

	QPainter painter(this);

	QPen penLine;
	penLine.setColor(Qt::yellow);
	penLine.setWidth(4);

	QPen penPoint;
	penPoint.setWidth(5);

	QLinearGradient gradient;
	painter.setBrush(gradient);
	QPointF pointPrevious;
	QPolygonF points;
	for(int n=0; n<m_variableDensity->steps; n++)
	{
		QPointF point(m_variableDensity->step[n].kr, m_variableDensity->step[n].scale);		QPointF pointCurrent = mapToWindow(point);
		points.append(pointCurrent);
		if(n)
		{
			painter.setPen(penLine);
			painter.drawLine(pointPrevious, pointCurrent);
		}
		pointPrevious = pointCurrent;
	}
	for(int n=0; n<m_variableDensity->steps; n++)
	{
		if(n==m_hoverIndex)
			penPoint.setColor(Qt::red);
		else
			penPoint.setColor(Qt::yellow);
		QPainterPath path;
		path.addEllipse(points.at(n), 4, 4);
		painter.setPen(penPoint);
		painter.drawPath(path);
	}
}

void VariableDensityPlot::mouseMoveEvent(QMouseEvent *event)
{
	QPointF cursorPoint = event->pos();
//	qWarning() << cursorPoint;

	m_hoverIndex = None;
	for(int n=0; n<m_variableDensity->steps; n++)
	{
		const struct VariableDensityStep* step = &m_variableDensity->step[n];
		QPointF variableDensityPoint = mapToWindow(QPointF(step->kr, step->scale));
		QPointF offset = cursorPoint - variableDensityPoint;
		float distance = offset.manhattanLength();

		if(distance<10)
		{
			m_hoverIndex = n;
//			qWarning() << distance;
			update();
			break;
		}
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

void VariableDensityPlot::addPoint(float kr, float scale)
{
	addLinearVariableDensityStep(m_variableDensity, kr, scale);
}


