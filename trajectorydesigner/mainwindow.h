#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

QT_FORWARD_DECLARE_CLASS(QSlider)
QT_FORWARD_DECLARE_CLASS(QSpinBox)
QT_FORWARD_DECLARE_CLASS(QDoubleSpinBox)
QT_FORWARD_DECLARE_CLASS(QwtPlot)
QT_FORWARD_DECLARE_CLASS(QwtPlotCurve)
QT_FORWARD_DECLARE_STRUCT(Trajectory)

namespace Ui {
class MainWindow;
}

QT_FORWARD_DECLARE_CLASS(Generator)

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	explicit MainWindow(QWidget *parent = 0);
	~MainWindow();
private slots:
	void updateFieldOfViewDisplay();
	void setFieldOfView(float fieldOfView, int axis);
	void setSpatialResolution(float spatialResolution, int axis);
	void setReadoutDuration(float duration);
	void setReadouts(int readouts);
	void updateTrajectoryPlot(Trajectory* trajectory);
	void updateGradientsPlot(Trajectory* trajectory);

private:
	Ui::MainWindow *ui;
	QSlider* m_fieldOfViewSlider[3];
	QSpinBox* m_fieldOfViewSpinBox[3];

	QSlider* m_spatialResolutionSlider[3];
	QDoubleSpinBox* m_spatialResolutionSpinBox[3];
	float m_spatialResolutionSliderScale = 10;
	float m_readoutDurationSliderScale = 10;

	QwtPlot* m_trajectoryPlot;
	QwtPlotCurve* m_trajectoryCurve;

	QwtPlot* m_gradientsPlot;
	QwtPlotCurve* m_gradientPlotCurve[3];

	Generator * m_generator;
};

#endif // MAINWINDOW_H
