#include "generator.h"

extern "C"
{
#include "spiral.h"
}

Generator::Generator(QObject *parent) : QObject(parent),
	m_trajectoryType(Spiral)
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

void Generator::setFieldOfView(float fieldOfView, int axis)
{
	m_fieldOfView[axis] = fieldOfView;
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
	switch(m_trajectoryType)
	{
		case Spiral:
			generateSpirals(NULL, m_fieldOfView[0], m_spatialResolution[0], m_readoutDuration, 4e-6, m_readouts, Archimedean, 0, m_fieldOfView[0], 4, 15000);
			break;
	}
}
