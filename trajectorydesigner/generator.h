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
	float m_fieldOfView[3];
	float m_spatialResolution[3];
	float m_readoutDuration;
	int m_readouts;
	Trajectory* m_trajectory;
};

Q_DECLARE_METATYPE(Generator::TrajectoryType)

#endif // GENERATOR_H
