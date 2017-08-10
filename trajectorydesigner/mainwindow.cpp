#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "generator.h"

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
	connect(ui->trajectoryComboBox, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged), [=](int index){
		Generator::TrajectoryType type = ui->trajectoryComboBox->currentData().value<Generator::TrajectoryType>();
		m_generator->setTrajectory(type);
	});

	updateFieldOfViewDisplay();
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::updateFieldOfViewDisplay()
{
	QWidget* layout;
	switch(m_generator->trajectoryType())
	{
		case Generator::Spiral:
//			layout = qobject_cast<QWidget*>(ui->fieldOfViewYLayout);
//			layout->hide();
			break;
	}
}
