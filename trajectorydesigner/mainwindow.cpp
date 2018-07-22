#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "generator.h"
#include "timeseriesplot.h"
#include "plot2d.h"
#include "variabledensitydesigner.h"
#include "phantomreconstruction.h"
extern "C"
{
#include "trajectory.h"
}

#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>

#include <QDebug>
#include <QCheckBox>

MainWindow::MainWindow(QWidget *parent) :
	QMainWindow(parent),
	ui(new Ui::MainWindow),
	m_trajectoryPlotXZ(new Plot2D(this)),
	m_slewRatePlot(new TimeSeriesPlot(3)),
	m_generator(new Generator),
	m_phantomReconstruction(new PhantomReconstruction())
{
  ui->setupUi(this);
  ui->trajectoryComboBox->addItem("Spiral", SPIRAL);
  ui->trajectoryComboBox->addItem("Radial", RADIAL);
  ui->trajectoryComboBox->addItem("Radial 3D", RADIAL3D);
  ui->trajectoryComboBox->addItem("Cones", CONES);
  ui->trajectoryComboBox->addItem("Rings", RINGS);

//	connect(ui->trajectoryComboBox, SIGNAL(currentIndexChanged(int)), m_generator, SLOT(setTrajectory(TrajectoryType)));
	connect(ui->trajectoryComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(updateFieldOfViewDisplay()));
//	connect(ui->trajectoryComboBox, SIGNAL(currentIndexChanged(int)), [=](int index){
//		Generator::TrajectoryType type = ui->trajectoryComboBox->currentData().value<Generator::TrajectoryType>();
//		m_generator->setTrajectory(type);
//	});
	connect(ui->trajectoryComboBox, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged), [=](int index){
    TrajectoryType type = static_cast<TrajectoryType>(ui->trajectoryComboBox->itemData(index).toInt());
		qWarning() << type;
		bool autoUpdate;
		switch(type) {
      case CONES:
				autoUpdate = false;
				break;
			default:
				autoUpdate = true;
				break;
		}
		setAutoUpdate(autoUpdate);
		m_generator->setTrajectory(type);
	});
	connect(ui->variableDensityPushButton, SIGNAL(clicked(bool)), this, SLOT(setVariableDensity(bool)));

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
		setFieldOfView(28, n);
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
		m_spatialResolutionSpinBox[n]->setMaximum(5);
		m_spatialResolutionSlider[n]->setMinimum(m_spatialResolutionSpinBox[n]->minimum()*m_spatialResolutionSliderScale);
		m_spatialResolutionSlider[n]->setMaximum(m_spatialResolutionSpinBox[n]->maximum()*m_spatialResolutionSliderScale);

		connect(m_spatialResolutionSlider[n], &QSlider::valueChanged, [=](int value){
			setSpatialResolution(value/m_spatialResolutionSliderScale, n);
		});
		connect(m_spatialResolutionSpinBox[n], &QSpinBox::editingFinished, [=]{
			setSpatialResolution(m_spatialResolutionSpinBox[n]->value(), n);
		});
		setSpatialResolution(2, n);
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
//  m_trajectoryPlot->setAxisAutoScale(QwtPlot::xBottom, false);
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

	m_slewRatePlot->setYLimits(-15000, 15000);

	ui->gridLayout->addWidget(m_trajectoryPlot, 0, 0);
	ui->gridLayout->addWidget(m_gradientsPlot, 0, 1);
	ui->gridLayout->addWidget(m_slewRatePlot, 1, 1);
	ui->gridLayout->addWidget(m_trajectoryPlotXZ, 1, 0);

	ui->tabWidget->setCurrentIndex(1);
	QPointer<QWidget> tab = ui->tabWidget->currentWidget();
  ui->phantomLayout->addWidget(m_phantomReconstruction);

	connect(m_generator, SIGNAL(updated(Trajectory*)), this, SLOT(updateTrajectoryPlot(Trajectory*)));
	connect(m_generator, SIGNAL(updated(Trajectory*)), this, SLOT(updateGradientsPlot(Trajectory*)));
	connect(m_generator, SIGNAL(updated(Trajectory*)), this, SLOT(updateSlewRatePlot(Trajectory*)));
  connect(m_generator, &Generator::updated, this, &MainWindow::updateReadoutIndex);
  connect(m_generator, SIGNAL(updated(Trajectory*)), m_phantomReconstruction.data(), SLOT(reconstruct(Trajectory*)));

	connect(m_generator, SIGNAL(readoutsChanged(int)), this, SLOT(setReadouts(int)));
	connect(ui->autoUpdateCheckBox, SIGNAL(toggled(bool)), m_generator, SLOT(setAutoUpdate(bool)));
	connect(ui->autoUpdateCheckBox, &QCheckBox::toggled, [=](bool status) {
		ui->updatePushButton->setEnabled(!status);
	});
	ui->autoUpdateCheckBox->setChecked(true);
	connect(ui->updatePushButton, SIGNAL(clicked()), m_generator, SLOT(update()));

	connect(m_generator, &Generator::updated, [=](Trajectory* trajectory) {
    setPlotReadouts(trajectory->numReadouts-1);
	});

  connect(ui->tabWidget, &QTabWidget::currentChanged, [=](int index)
  {
    if(index==1)
    {
      m_phantomReconstruction->setEnabled(true);
      m_phantomReconstruction->reconstruct(m_generator->trajectory());
    }
    else
      m_phantomReconstruction->setEnabled(false);
  });
	ui->tabWidget->setCurrentIndex(0);

  connect(ui->readoutSlider, &QSlider::valueChanged, this, &MainWindow::setReadout);
  connect(ui->readoutSpinBox, &QSpinBox::editingFinished, [=]
    {
      setReadout(ui->readoutSpinBox->value());
    }
  );

  m_generator->update();
}

MainWindow::~MainWindow()
{
	delete ui;
	delete m_trajectoryPlotXZ;
}

void MainWindow::updateFieldOfViewDisplay()
{
	switch(m_generator->trajectoryType())
	{
    case SPIRAL:
    case RINGS:
//			layout = qobject_cast<QWidget*>(ui->fieldOfViewYLayout);
//			layout->hide();
			break;
    case RADIAL:
			break;
    case CONES:
			break;
    case RADIAL3D:
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
	QVector<float> spatialResolutions = m_generator->spatialResolution();
	float minSpatialResolution = spatialResolutions[0];
	float exent = 10/minSpatialResolution;
	if(m_trajectoryPlot)
	{
		m_trajectoryPlot->setAxisScale(QwtPlot::xBottom, -exent/2, exent/2);
		m_trajectoryPlot->setAxisScale(QwtPlot::yLeft, -exent/2, exent/2);
	}
	if(m_trajectoryPlotXZ)
	{
		m_trajectoryPlotXZ->setAxisScale(QwtPlot::xBottom, -exent/2, exent/2);
		m_trajectoryPlotXZ->setAxisScale(QwtPlot::yLeft, -exent/2, exent/2);
	}
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

void MainWindow::setPlotReadouts(int readouts)
{
  ui->readoutsSlider->setMaximum(readouts);
  ui->readoutSpinBox->setMaximum(readouts);
}

void MainWindow::setReadout(int readout)
{
  ui->readoutsSlider->setValue(readout);
  ui->readoutSpinBox->setValue(readout);
	bool changed = m_readout != readout;
	m_readout = readout;

	Trajectory* trajectory = m_generator->trajectory();
	if(changed && trajectory)
	{
		updateTrajectoryPlot(trajectory);
		updateGradientsPlot(trajectory);
		updateSlewRatePlot(trajectory);
	}
}

void MainWindow::setAutoUpdate(bool status)
{
	m_generator->setAutoUpdate(status);
	ui->autoUpdateCheckBox->setChecked(status);
}

void MainWindow::setVariableDensity(bool status)
{
	m_generator->setVariableDensity(status);
	if(status)
	{
		if(!m_variableDensityDesigner)
		{
      m_variableDensityDesigner = new VariableDensityDesigner(m_generator->variableDensity());
			connect(m_variableDensityDesigner, &VariableDensityDesigner::updated, [=]()
			{
				if(m_generator->autoUpdate())
					m_generator->update();
      });
		}
		m_variableDensityDesigner->show();
	}
	else
	{
		m_variableDensityDesigner->hide();
	}
}

void MainWindow::updateTrajectoryPlot(Trajectory *trajectory)
{
	QVector<QPointF> coordinatesXY;
	QVector<QPointF> coordinatesXZ;
	float kSpaceCoordinates[3];
  for(int n=0; n<trajectory->numReadoutPoints; n++)
	{
		trajectoryCoordinates(n, m_readout, trajectory, kSpaceCoordinates, NULL);
		coordinatesXY.append(QPointF(kSpaceCoordinates[0], kSpaceCoordinates[1]));
    if(trajectory->numDimensions>2)
      coordinatesXZ.append(QPointF(kSpaceCoordinates[0], kSpaceCoordinates[2]));
	}

	m_trajectoryCurve->setSamples(coordinatesXY);
	m_trajectoryPlotXZ->setSamples(coordinatesXZ);

//	m_trajectoryPlot->resize(300,300);
	m_trajectoryPlot->replot();
	m_trajectoryPlotXZ->replot();
}

void MainWindow::updateGradientsPlot(Trajectory *trajectory)
{
	QVector<QPointF> coordinates;
  for(int d=0; d<trajectory->numDimensions; d++)
	{
		coordinates.clear();
		float* waveform = trajectoryGradientWaveform(trajectory, m_readout, d);
    for(int n=0; n<trajectory->numWaveformPoints; n++)
		{
			float t = n*trajectory->samplingInterval*1e3;
			coordinates.append(QPointF(t,waveform[n]));
		}
		m_gradientPlotCurve[d]->setSamples(coordinates);
	}

	m_gradientsPlot->replot();
}

void MainWindow::updateSlewRatePlot(Trajectory *trajectory)
{
  m_slewRatePlot->setTime(trajectory->samplingInterval, trajectory->numWaveformPoints-1);
  for(int d=0; d<trajectory->numDimensions; d++)
	{
		float* gradientWaveform = trajectoryGradientWaveform(trajectory, m_readout, d);
		QVector<double> slew;
    for(int n=0; n<trajectory->numWaveformPoints-1; n++)
			slew.append((gradientWaveform[n+1]-gradientWaveform[n])/trajectory->samplingInterval);
		m_slewRatePlot->setSeriesData(slew, d);
	}
  m_slewRatePlot->replot();
}

void MainWindow::updateReadoutIndex(Trajectory *trajectory)
{
  int maxIndex = trajectory->storage==STORE_ALL ? trajectory->numReadouts : trajectory->numBases;
  maxIndex--;
  ui->readoutSlider->setMaximum(maxIndex);
  ui->readoutSpinBox->setMaximum(maxIndex);
}

