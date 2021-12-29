#include "trajectorygenerator.h"

#include <QCoreApplication>
#include <QCommandLineParser>
#include <QDebug>

#include <iostream>


int main(int argc, char *argv[])
{
  QCoreApplication application(argc, argv);

  QCommandLineParser parser;
  parser.addHelpOption();
  parser.addOptions({{"trajectory", "The type of trajectory to generate", "trajectory"},
                     {"readdur", "Readout duration", "readoutDuration"},
                     {"fov", "Field of view", "fieldOfView"},
                     {"res", "Spatial resolution", "spatialResolution"},
                     {"comp", "Cones compensation", "conesCompensation"},
                     {"bases", "Number of basis waveforms", "numBases"},
                     {"to", "File path to save the trajectory to", "outputFilePath"}});
  parser.process(application);

  if(parser.positionalArguments().empty())
  {
    qCritical("No arguments specified");
    return 1;
  }

  TrajectoryType trajectoryType;

  const QString trajectoryTypeRequested = parser.positionalArguments()[0];
  if(!trajectoryTypeRequested.compare("cones", Qt::CaseInsensitive))
     trajectoryType = CONES;
  else
    qCritical("Invalid trajectory type %s", qPrintable(trajectoryTypeRequested));

  TrajectoryGenerator generator;
  generator.setTrajectoryType(trajectoryType);
  if(parser.isSet("readdur"))
    generator.setReadoutDuration(parser.value("readdur").toFloat());

  if(parser.isSet("fov"))
  {
    const QStringList fieldOfView = parser.values("fov");
    if (trajectoryType == CONES)
    {
      const float fieldOfViewXY = fieldOfView[0].toFloat();
      const float fieldOfViewZ = fieldOfView[1].toFloat();
      generator.setFieldOfView({fieldOfViewXY, fieldOfViewXY, fieldOfViewZ});
    }
  }

  if(parser.isSet("res"))
  {
    const QStringList spatialResolution = parser.values("res");
    if (trajectoryType == CONES)
    {
      const float spatialResolutionXY = spatialResolution[0].toFloat();
      const float spatialResolutionZ = spatialResolution[1].toFloat();
      generator.setSpatialResolution({spatialResolutionXY, spatialResolutionXY, spatialResolutionZ});
    }
  }

  if(parser.isSet("bases"))
    generator.setNumBases(parser.value("bases").toInt());

  if(parser.isSet("comp") && trajectoryType == CONES)
    generator.setCompensation(parser.value("comp").toStdString());

  std::cout << generator.info();
  generator.generate();

  if(parser.isSet("to"))
  {
    const QString filePath = parser.value("to");
    qInfo() << "Saving trajectory to " << filePath;
    generator.save(filePath.toStdString());
  }

  return 0;
}
