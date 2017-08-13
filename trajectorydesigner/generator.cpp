#include "generator.h"

extern "C"
{
#include "spiral.h"
#include "radial.h"
}

Generator::Generator(QObject *parent) : QObject(parent)
{

}

Generator::TrajectoryType Generator::trajectoryType()
{
	return m_trajectoryType;
}

void Generator::setTrajectory(TrajectoryType type)
{
	bool changed = m_trajectoryType != type;
	m_trajectoryType = type;
	if(changed)
		update();
}

void Generator::setFieldOfView(float fieldOfView, int axis)
{
	bool changed = m_fieldOfView[axis] != fieldOfView;
	m_fieldOfView[axis] = fieldOfView;

	if(changed)
		update();
}

void Generator::setSpatialResolution(float spatialResolution, int axis)
{
	bool changed = m_spatialResolution[axis] != spatialResolution;
	m_spatialResolution[axis] = spatialResolution;

	if(changed)
		update();
}

void Generator::setReadoutDuration(float readoutDuration)
{
	bool changed = m_readoutDuration != readoutDuration;
	m_readoutDuration = readoutDuration;

	if(changed)
		update();
}

void Generator::update()
{
	int previousReadouts = m_trajectory ? m_trajectory->readouts : 0;

	switch(m_trajectoryType)
	{
		case Spiral:
			m_trajectory = generateSpirals(NULL, m_fieldOfView[0], m_spatialResolution[0], m_readoutDuration, 4e-6, m_readouts, Archimedean, 0, m_fieldOfView[0], 4, 15000);
			break;
		case Radial:
		m_trajectory = generateRadial2D(m_fieldOfView[0], m_fieldOfView[1], EllipticalShape, m_spatialResolution[0], m_spatialResolution[1], EllipticalShape, 1, 1, 4, 15000, 4e-6);
			break;
	}

	emit updated(m_trajectory);

	if(previousReadouts!=m_trajectory->readouts)
		emit readoutsChanged(m_trajectory->readouts);
}
