#ifndef VARIABLEDENSITYDESIGNER_H
#define VARIABLEDENSITYDESIGNER_H

#include <QWidget>
#include <QPointer>

namespace Ui {
class VariableDensityDesigner;
}

QT_FORWARD_DECLARE_CLASS(VariableDensityPlot)

class VariableDensityDesigner : public QWidget
{
	Q_OBJECT

public:
	explicit VariableDensityDesigner(QWidget *parent = 0);
	~VariableDensityDesigner();

private:
	Ui::VariableDensityDesigner *ui;
	QPointer<VariableDensityPlot> m_variableDensityPlot;
};

#endif // VARIABLEDENSITYDESIGNER_H
