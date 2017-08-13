#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "generator.h"
extern "C"
{
#include "trajectory.h"
}

#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>

MainWindow::MainWindow(QWidget *parent) :
	QMainWindow(parent),
	ui(new Ui::MainWindow),
	m_generator(new Generator)
{
	ui->setupUi(this);
	ui->trajectoryComboBox->addItem("Spiral", Generator::Spiral);

//	connect(ui->trajectoryComboBox, SIGNAL(currentIndexChanged(int)), m_generator, SLOT(setTrajectory(TrajectoryType)));
	connect(ui->trajectoryComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(updateFieldOfViewDisplay()));
//	connect(ui->trajectoryComboBox, SIGNAL(currentIndexChanged(int)), [=](int index){
//		Generator::TrajectoryType type = ui->trajectoryComboBox->currentData().value<Generator::TrajectoryType>();
//		m_generator->setTrajectory(type);
//	});
	connect(ui->trajectoryComboBox, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged), [=](){
		Generator::TrajectoryType type = ui->trajectoryComboBox->currentData().value<Generator::TrajectoryType>();
		m_generator->setTrajectory(type);
	});

	updateFieldOfViewDisplay();

	m_fieldOfViewSlider[0] = ui->fieldOfViewXSlider;
	m_fieldOfViewSlider[1] = ui->fieldOfViewYSlider;
	m_fieldOfViewSlider[2] = ui->fieldOfViewZSlider;
	m_fieldOfViewSpinBox[0] = ui->fieldOfViewXSpinBox;
	m_fieldOfViewSpinBox[1] = ui->fieldOfViewYSpinBox;
	m_fieldOfViewSpinBox[2] = ui->fieldOfViewZSpinBox;

	for (int n=0; n<3; n++)
	{
		m_fieldOfViewSlider[n]->setMinimum(5);
		m_fieldOfViewSlider[n]->setMaximum(50);
		m_fieldOfViewSpinBox[n]->setMinimum(m_fieldOfViewSlider[n]->minimum());
		m_fieldOfViewSpinBox[n]->setMaximum(m_fieldOfViewSlider[n]->maximum());
		connect(m_fieldOfViewSlider[n], &QSlider::valueChanged, [=](int value){
			setFieldOfView(value, n);
		});
		connect(m_fieldOfViewSpinBox[n], &QSpinBox::editingFinished, [=]{
			setFieldOfView(m_fieldOfViewSpinBox[n]->value(), n);
		});
		m_fieldOfViewSlider[n]->setValue(28);
	}

	m_spatialResolutionSlider[0] = ui->spatialResolutionXSlider;
	m_spatialResolutionSlider[1] = ui->spatialResolutionYSlider;
	m_spatialResolutionSlider[2] = ui->spatialResolutionZSlider;
	m_spatialResolutionSpinBox[0] = ui->spatialResolutionXDoubleSpinBox;
	m_spatialResolutionSpinBox[1] = ui->spatialResolutionYDoubleSpinBox;
	m_spatialResolutionSpinBox[2] = ui->spatialResolutionZDoubleSpinBox;

	for (int n=0; n<3; n++)
	{
		m_spatialResolutionSpinBox[n]->setMinimum(.3);
		m_spatialResolutionSpinBox[n]->setMaximum(12);
		m_spatialResolutionSlider[n]->setMinimum(m_spatialResolutionSpinBox[n]->minimum()*m_spatialResolutionSliderScale);
		m_spatialResolutionSlider[n]->setMaximum(m_spatialResolutionSpinBox[n]->maximum()*m_spatialResolutionSliderScale);

		connect(m_spatialResolutionSlider[n], &QSlider::valueChanged, [=](int value){
			setSpatialResolution(value/m_spatialResolutionSliderScale, n);
		});
		connect(m_spatialResolutionSpinBox[n], &QSpinBox::editingFinished, [=]{
			setSpatialResolution(m_spatialResolutionSpinBox[n]->value(), n);
		});
		m_spatialResolutionSpinBox[n]->setValue(2);
	}

	ui->readoutDurationDoubleSpinBox->setMinimum(.128);
	ui->readoutDurationDoubleSpinBox->setMaximum(20);
	ui->readoutDurationSlider->setMinimum(ui->readoutDurationDoubleSpinBox->minimum()*m_readoutDurationSliderScale);
	ui->readoutDurationSlider->setMaximum(ui->readoutDurationDoubleSpinBox->maximum()*m_readoutDurationSliderScale);
	connect(ui->readoutDurationSlider, &QSlider::valueChanged, [=](int value){
		setReadoutDuration(value/m_readoutDurationSliderScale);
	});
	connect(ui->readoutDurationDoubleSpinBox, &QSpinBox::editingFinished, [=]{
		setReadoutDuration(ui->readoutDurationDoubleSpinBox->value());
	});
	setReadoutDuration(8);

	ui->readoutsSlider->setEnabled(false);
	ui->readoutsSpinBox->setEnabled(false);
//	ui->readoutsSlider->setMinimum(1);
//	ui->readoutsSpinBox->setMinimum(ui->readoutsSlider->minimum());
//	connect(ui->readoutsSlider, &QSlider::valueChanged, [=](int value){
//		setFieldOfView(value, n);
//	});

	m_trajectoryPlot = new QwtPlot(QwtText("Trajectory Plot"), parent);
	m_trajectoryPlot->setCanvasBackground(Qt::black);
	QwtPlotGrid *trajectoryPlotGrid = new QwtPlotGrid;
	trajectoryPlotGrid->attach(m_trajectoryPlot);
	trajectoryPlotGrid->setMajorPen(Qt::gray, 0, Qt::DotLine);
	m_trajectoryCurve = new QwtPlotCurve;
	m_trajectoryCurve->setPen(Qt::green);
	m_trajectoryCurve->attach(m_trajectoryPlot);

	m_gradientsPlot = new QwtPlot(QwtText("Gradients"), parent);
	m_gradientsPlot->setCanvasBackground(Qt::black);
	m_gradientsPlot->setAxisScale(QwtPlot::yLeft, -4, 4);
	m_gradientsPlot->setAxisScale(QwtPlot::xBottom, 0, 20);
	for(int n=0; n<3; n++)
	{
		m_gradientPlotCurve[n] = new QwtPlotCurve;
		m_gradientPlotCurve[n]->attach(m_gradientsPlot);
	}
	m_gradientPlotCurve[0]->setPen(Qt::blue);
	m_gradientPlotCurve[1]->setPen(Qt::green);

	ui->gridLayout->addWidget(m_trajectoryPlot, 0, 0);
	ui->gridLayout->addWidget(m_gradientsPlot, 0, 1);

	connect(m_generator, SIGNAL(updated(Trajectory*)), this, SLOT(updateTrajectoryPlot(Trajectory*)));
	connect(m_generator, SIGNAL(updated(Trajectory*)), this, SLOT(updateGradientsPlot(Trajectory*)));

	connect(m_generator, SIGNAL(readoutsChanged(int)), this, SLOT(setReadouts(int)));
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::updateFieldOfViewDisplay()
{
	//QWidget* layout;
	switch(m_generator->trajectoryType())
	{
		case Generator::Spiral:
//			layout = qobject_cast<QWidget*>(ui->fieldOfViewYLayout);
//			layout->hide();
			break;
	}
}

void MainWindow::setFieldOfView(float fieldOfView, int axis)
{
	m_fieldOfViewSlider[axis]->setValue(fieldOfView);
	m_fieldOfViewSpinBox[axis]->setValue(fieldOfView);
	m_generator->setFieldOfView(fieldOfView, axis);
}

void MainWindow::setSpatialResolution(float spatialResolution, int axis)
{
	m_spatialResolutionSlider[axis]->setValue(m_spatialResolutionSliderScale*spatialResolution);
	m_spatialResolutionSpinBox[axis]->setValue(spatialResolution);
	m_generator->setSpatialResolution(spatialResolution, axis);
}

void MainWindow::setReadoutDuration(float duration)
{
	ui->readoutDurationSlider->setValue(duration*m_readoutDurationSliderScale);
	ui->readoutDurationDoubleSpinBox->setValue(duration);
	m_generator->setReadoutDuration(duration*1e-3);
}

void MainWindow::setReadouts(int readouts)
{
	ui->readoutsSlider->setValue(readouts);
	ui->readoutsSpinBox->setValue(readouts);
}

void MainWindow::updateTrajectoryPlot(Trajectory *trajectory)
{
	QVector<QPointF> coordinates;
	float kSpaceCoordinates[2];
	for(int n=0; n<trajectory->readoutPoints; n++)
	{
		trajectoryCoordinates(n, 0, trajectory, kSpaceCoordinates, NULL);
		coordinates.append(QPointF(kSpaceCoordinates[0], kSpaceCoordinates[1]));
	}

	m_trajectoryCurve->setSamples(coordinates);

//	m_trajectoryPlot->resize(300,300);
	m_trajectoryPlot->replot();
}

void MainWindow::updateGradientsPlot(Trajectory *trajectory)
{
	QVector<QPointF> coordinates;
	for(int d=0; d<trajectory->dimensions; d++)
	{
		coordinates.clear();
		float* waveform = trajectoryGradientWaveform(trajectory, 0, d);
		for(int n=0; n<trajectory->waveformPoints; n++)
		{
			float t = n*trajectory->samplingInterval*1e3;
			coordinates.append(QPointF(t,waveform[n]));
		}
		m_gradientPlotCurve[d]->setSamples(coordinates);
	}

	m_gradientsPlot->replot();
}
