#ifndef VARIABLEDENSITYPLOT_H
#define VARIABLEDENSITYPLOT_H

#include <QWidget>

QT_FORWARD_DECLARE_STRUCT(VariableDensity)

class VariableDensityPlot : public QWidget
{
	Q_OBJECT
public:
	explicit VariableDensityPlot(struct VariableDensity* variableDensity, QWidget *parent = 0);

signals:
	void updated();
public slots:
protected:
	void paintEvent(QPaintEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mousePressEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
private:
	enum Indices{None=-1};

	QPointF mapToWindow(QPointF point);
	QPointF mapToDensity(QPointF point);
	void drawLine(QPointF point1, QPointF point2, const QPen& pen);
	void addPoint(float kr, float scale);
	float offsetFactor();

	QVector<QPointF> m_points = {{0,1}, {1,1}};
	VariableDensity* m_variableDensity = NULL;
	int m_hoverIndex = None;
	bool m_clicked = false;

	const float m_scale = .9;
};

#endif // VARIABLEDENSITYPLOT_H
