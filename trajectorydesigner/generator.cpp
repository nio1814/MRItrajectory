#include "generator.h"

extern "C"
{
#include "spiral.h"
#include "radial.h"
#include "cones.h"
#include "variabledensity.h"
}

Generator::Generator(QObject *parent) : QObject(parent)
{

}


void Generator::setTrajectory(TrajectoryType type)
{
	bool changed = m_trajectoryType != type;
	m_trajectoryType = type;
	TrajectoryGenerator::setTrajectoryType(type);

	if(changed && m_autoUpdate)
		update();
}

void Generator::setFieldOfView(float fieldOfView, int axis)
{
	bool changed = m_fieldOfView[axis] != fieldOfView;
	m_fieldOfView[axis] = fieldOfView;

	if(changed && m_autoUpdate)
		update();
}

QVector<float> Generator::fieldOfView()
{
	QVector<float> fov(3);
	memcpy(fov.data(), m_fieldOfView, 3*sizeof(float));

	return fov;
}

void Generator::setSpatialResolution(float spatialResolution, int axis)
{
	bool changed = m_spatialResolution[axis] != spatialResolution;
	m_spatialResolution[axis] = spatialResolution;

	if(changed && m_autoUpdate)
		update();
}

QVector<float> Generator::spatialResolution()
{
	QVector<float> res(3);
	memcpy(res.data(), m_spatialResolution, 3*sizeof(float));

	return res;
}

void Generator::setReadoutDuration(float readoutDuration)
{
	bool changed = m_readoutDuration != readoutDuration;
	m_readoutDuration = readoutDuration;

	if(changed && m_autoUpdate)
		update();
}

void Generator::setAutoUpdate(bool status)
{
	m_autoUpdate = status;
}

bool Generator::autoUpdate()
{
	return m_autoUpdate;
}

void Generator::setVariableDensity(bool status)
{
	if(status)
	{
		if(!m_variableDensity)
		{
			m_variableDensity = newVariableDensity();
			setToFullySampled(m_variableDensity);
		}
	}
	else
		deleteVariableDensity(&m_variableDensity);
}

VariableDensity *Generator::variableDensity()
{
	return m_variableDensity;
}

void Generator::update()
{
  int previousReadouts = m_trajectory ? m_trajectory->numReadouts : 0;
	generate();

	emit updated(m_trajectory);

  if(previousReadouts!=m_trajectory->numReadouts)
    emit readoutsChanged(m_trajectory->numReadouts);
}
