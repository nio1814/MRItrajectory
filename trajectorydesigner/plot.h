#ifndef PLOT_H
#define PLOT_H

#include <qwt_plot.h>

QT_FORWARD_DECLARE_CLASS(QwtPlotCurve)

class Plot : public QwtPlot
{
	Q_OBJECT
public:
	explicit Plot(int curves=1, QWidget *parent = 0);
	void setYLimits(float yMin, float yMax);
signals:

public slots:
protected:
	QVector<QSharedPointer<QwtPlotCurve> > m_plotCurve;
private:
	QColor m_colorOrder[3] = {Qt::yellow, QColor(255,0,255), QColor(0,255,255)};
};

#endif // PLOT_H
