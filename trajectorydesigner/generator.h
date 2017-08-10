#ifndef GENERATOR_H
#define GENERATOR_H

#include <QObject>

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
	void setFieldOfView(float value, int axis);
private:
	TrajectoryType m_trajectoryType;
	float m_fieldOfView[3];
};

Q_DECLARE_METATYPE(Generator::TrajectoryType)

#endif // GENERATOR_H
