from argparse import ArgumentParser
from sys import argv

from cmri.cmri import TrajectoryGenerator, TrajectoryType


if __name__ == "__main__":
  parser = ArgumentParser()
  parser.add_argument('type', type=str)
  parser.add_argument('--readdur', type=float, dest='readout_duration')
  parser.add_argument('--fov', type=float, nargs='+', dest='field_of_view')
  parser.add_argument('--res', type=float, nargs='+', dest='spatial_resolution')
  parser.add_argument('--bases', type=int)
  parser.add_argument('--to', type=str, dest='file_path_output')
  arguments = parser.parse_args(argv[1:])

  trajectory_type = arguments.type.lower()
  if trajectory_type == 'cones':
    trajectory_type = TrajectoryType.CONES
  else:
    raise ValueError(f'Unrecognized trajectory type {trajectory_type}')

  generator = TrajectoryGenerator()
  generator.setTrajectoryType(trajectory_type)
  if arguments.readout_duration:
    generator.setReadoutDuration(arguments.readout_duration)
  if arguments.field_of_view:
    generator.setFieldOfView(arguments.field_of_view)
  if arguments.spatial_resolution:
    generator.setSpatialResolution(arguments.spatial_resolution)
  if arguments.bases:
    generator.setNumBases(arguments.bases)
  
  generator.generate()

  if arguments.file_path_output:
    print(f'Saving trajectory to {arguments.file_path_output}')
    generator.save(arguments.file_path_output)