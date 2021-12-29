#include "variabledensitydesigner.h"
#include "ui_variabledensitydesigner.h"

#include "variabledensityplot.h"

VariableDensityDesigner::VariableDensityDesigner(VariableDensity *variableDensity, QWidget *parent) :
	QWidget(parent),
	ui(new Ui::VariableDensityDesigner),
	m_variableDensityPlot(new VariableDensityPlot(variableDensity, this))
{
	ui->setupUi(this);
	ui->verticalLayout->addWidget(m_variableDensityPlot);

	connect(m_variableDensityPlot, SIGNAL(updated()), this, SIGNAL(updated()));
}

VariableDensityDesigner::~VariableDensityDesigner()
{
	delete ui;
}
