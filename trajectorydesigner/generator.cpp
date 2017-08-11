#include "generator.h"

#include "spiral.h"

Generator::Generator(QObject *parent) : QObject(parent),
	m_trajectoryType(Spiral),
	m_trajectory(m_trajectory)
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

void Generator::update()
{
	switch(m_trajectoryType)
	{
		case Spiral:
			break;
	}
}
