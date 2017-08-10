#include "generator.h"

Generator::Generator(QObject *parent) : QObject(parent), m_trajectoryType(Spiral)
{

}

Generator::TrajectoryType Generator::trajectoryType()
{
	return m_trajectoryType;
}

void Generator::setTrajectory(TrajectoryType type)
{
	m_trajectoryType = type;
}

void Generator::setFieldOfView(float value, int axis)
{
	m_fieldOfView[axis] = value;
}
