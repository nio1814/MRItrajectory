#ifndef TIMESERIESPLOT_H
#define TIMESERIESPLOT_H

#include "plot.h"

class TimeSeriesPlot : public Plot
{
public:
	TimeSeriesPlot(int numSeries, QWidget *parent=NULL);

	void setTime(double timeInterval, int points);
	void setSeriesData(float *data, int points, int index);
	void setSeriesData(QVector<double> data, int index);
private:
	double m_timeInterval = 1;
	QVector<double> m_time;
};

#endif // TIMESERIESPLOT_H
