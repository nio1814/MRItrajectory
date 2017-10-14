#include "variabledensitydesigner.h"
#include "ui_variabledensitydesigner.h"

#include "variabledensityplot.h"

VariableDensityDesigner::VariableDensityDesigner(QWidget *parent) :
	QWidget(parent),
	ui(new Ui::VariableDensityDesigner),
	m_variableDensityPlot(new VariableDensityPlot(this))
{
	ui->setupUi(this);
	ui->verticalLayout->addWidget(m_variableDensityPlot);
}

VariableDensityDesigner::~VariableDensityDesigner()
{
	delete ui;
}
