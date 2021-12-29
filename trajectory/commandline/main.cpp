#include <QCoreApplication>
#include <QCommandLineParser>
#include <QDebug>

#include "trajectorygenerator.h"


int main(int argc, char *argv[])
{
  QCoreApplication application(argc, argv);

  QCommandLineParser parser;
  parser.addHelpOption();
  parser.addOptions({{"trajectory", "The type of trajectory to generate"},
                     {"readdur", "Readout duration", "readoutDuration"},
                     {"fov", "Field of view", "fieldOfView"}});
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
  if(parser.isSet("readoutDuration"))
    generator.setReadoutDuration(parser.value("readoutDuration").toFloat());

  if(parser.isSet("fieldOfView"))
  {
    const QStringList fieldsOfView = parser.values("fieldOfView");
    qDebug() << fieldsOfView;
    if (trajectoryType == CONES)
    {
      const float fieldOfViewXY = fieldsOfView[0].toFloat();
      const float fieldOfViewZ = fieldsOfView[1].toFloat();
      generator.setFieldOfView({fieldOfViewXY, fieldOfViewXY, fieldOfViewZ});
    }
  }

  generator.generate();

  return 0;
}
