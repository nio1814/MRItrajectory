#ifndef GENERATOR_H
#define GENERATOR_H

#include <QObject>

QT_FORWARD_DECLARE_STRUCT(Trajectory)
QT_FORWARD_DECLARE_STRUCT(Cones)
QT_FORWARD_DECLARE_STRUCT(VariableDensity)

class Generator : public QObject
{
	Q_OBJECT
public:
	enum TrajectoryType{Spiral, Radial, Cones3D};
	explicit Generator(QObject *parent = 0);

	TrajectoryType trajectoryType();

signals:
	void updated(Trajectory *trajectory);
	void readoutsChanged(int readouts);
	void autoUpdate(bool status);

public slots:
	void setTrajectory(TrajectoryType type);
	Trajectory* trajectory();
	void setFieldOfView(float fieldOfView, int axis);
	void setSpatialResolution(float spatialResolution, int axis);
	void setReadoutDuration(float readoutDuration);
	void setAutoUpdate(bool status);
	bool autoUpdate();
	void setVariableDensity(bool status);
	VariableDensity *variableDensity();
	void update();
private:
	bool m_autoUpdate = true;
	TrajectoryType m_trajectoryType = Spiral;
	float m_fieldOfView[3] = {28,28,28};
	float m_spatialResolution[3] = {2,2,2};
	float m_readoutDuration = 8e-3;
	float m_gradientLimit = 4;
	float m_slewRateLimit = 15000;
	float m_samplingInterval = 4e-6;
	int m_fullProjection = 1;
	int m_readouts;
	VariableDensity* m_variableDensity = NULL;
	Trajectory* m_trajectory = NULL;
	Cones* m_cones = NULL;
};

Q_DECLARE_METATYPE(Generator::TrajectoryType)

#endif // GENERATOR_H
