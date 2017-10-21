#ifndef VARIABLEDENSITYDESIGNER_H
#define VARIABLEDENSITYDESIGNER_H

#include <QWidget>
#include <QPointer>

namespace Ui {
class VariableDensityDesigner;
}

QT_FORWARD_DECLARE_CLASS(VariableDensity)
QT_FORWARD_DECLARE_CLASS(VariableDensityPlot)

class VariableDensityDesigner : public QWidget
{
	Q_OBJECT

public:
	explicit VariableDensityDesigner(VariableDensity* variableDensity, QWidget *parent = 0);
	~VariableDensityDesigner();
signals:
	void updated();
private:
	Ui::VariableDensityDesigner *ui;
	QPointer<VariableDensityPlot> m_variableDensityPlot;
};

#endif // VARIABLEDENSITYDESIGNER_H
