#ifndef VARIABLEDENSITYPLOT_H
#define VARIABLEDENSITYPLOT_H

#include <QWidget>

class VariableDensityPlot : public QWidget
{
	Q_OBJECT
public:
	explicit VariableDensityPlot(QWidget *parent = 0);

signals:

public slots:
protected:
	void paintEvent(QPaintEvent *event);
private:
	QPointF mapToWindow(QPointF point);
	void drawLine(QPointF point1, QPointF point2, const QPen& pen);

	QVector<QPointF> m_points = {{0,1}, {1,1}};
};

#endif // VARIABLEDENSITYPLOT_H
