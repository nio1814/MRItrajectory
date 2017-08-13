#include "timeseriesplot.h"

#include <qwt_plot_curve.h>

TimeSeriesPlot::TimeSeriesPlot(int numSeries, QWidget* parent) :
	Plot(numSeries, parent)
{

}

void TimeSeriesPlot::setTime(double timeInterval, int points)
{
	m_timeInterval = timeInterval;
	m_time.clear();
	for(int n=0; n<points; n++)
		m_time.append(n*m_timeInterval);
}

void TimeSeriesPlot::setSeriesData(float *data, int points, int index)
{
	QVector<double> y(points);
	for(int n=0; n<points; n++)
		y[n] = data[n];
	m_plotCurve[index]->setSamples(m_time, y);
}

void TimeSeriesPlot::setSeriesData(QVector<double> data, int index)
{
	m_plotCurve[index]->setSamples(m_time, data);
}

