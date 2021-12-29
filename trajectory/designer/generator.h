#ifndef GENERATOR_H
#define GENERATOR_H

#include "trajectorygenerator.h"

#include <QObject>
#include <QVector>

QT_FORWARD_DECLARE_STRUCT(Trajectory)
QT_FORWARD_DECLARE_STRUCT(Cones)
QT_FORWARD_DECLARE_STRUCT(VariableDensity)

class Generator : public QObject, public TrajectoryGenerator
{
	Q_OBJECT
public:
	explicit Generator(QObject *parent = 0);

signals:
	void updated(Trajectory *trajectory);
	void readoutsChanged(int readouts);
	void autoUpdate(bool status);

public slots:
	void setTrajectory(TrajectoryType type);
	void setFieldOfView(float fieldOfView, int axis);
	QVector<float> fieldOfView();
	void setSpatialResolution(float spatialResolution, int axis);
	QVector<float> spatialResolution();
	void setReadoutDuration(float readoutDuration);
	void setAutoUpdate(bool status);
	bool autoUpdate();
	void setVariableDensity(bool status);
	VariableDensity *variableDensity();
	void update();
private:
	bool m_autoUpdate = true;
};

#endif // GENERATOR_H
