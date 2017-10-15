#ifndef VARIABLEDENSITYPLOT_H
#define VARIABLEDENSITYPLOT_H

#include <QWidget>

QT_FORWARD_DECLARE_STRUCT(VariableDensity)

class VariableDensityPlot : public QWidget
{
	Q_OBJECT
public:
	explicit VariableDensityPlot(QWidget *parent = 0);

signals:

public slots:
protected:
	void paintEvent(QPaintEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
private:
	enum Indices{None=-1};

	QPointF mapToWindow(QPointF point);
	void drawLine(QPointF point1, QPointF point2, const QPen& pen);
	void addPoint(float kr, float scale);

	QVector<QPointF> m_points = {{0,1}, {1,1}};
	VariableDensity* m_variableDensity = NULL;
	int m_hoverIndex = None;
};

#endif // VARIABLEDENSITYPLOT_H
