#ifndef GENERATOR_H
#define GENERATOR_H

#include <QObject>

QT_FORWARD_DECLARE_STRUCT(Trajectory)

class Generator : public QObject
{
	Q_OBJECT
public:
	enum TrajectoryType{Spiral};
	explicit Generator(QObject *parent = 0);

	TrajectoryType trajectoryType();

signals:

public slots:
	void setTrajectory(TrajectoryType type);
	void setFieldOfView(float fieldOfView, int axis);
	void setSpatialResolution(float spatialResolution, int axis);
	void setReadoutDuration(float readoutDuration);
private slots:
	void update();
private:
	TrajectoryType m_trajectoryType;
	float m_fieldOfView[3] = {28,28,28};
	float m_spatialResolution[3] = {2,2,2};
	float m_readoutDuration = 8e-3;
	int m_readouts;
	Trajectory* m_trajectory = NULL;
};

Q_DECLARE_METATYPE(Generator::TrajectoryType)

#endif // GENERATOR_H
