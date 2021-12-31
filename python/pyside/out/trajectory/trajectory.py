from enum import auto, Enum

import numpy as np

def read_array(file, data_type) -> np.ndarray:
  """Read an array stored in a binary file.

  Parameters
  ----------
  file : 
      The data file.
  data_type : 
      The data type of the stored array.
  """
  num_points = np.fromfile(file, 'i8', 1)[0]
  return np.fromfile(file, data_type, num_points)

class TrajectoryType(Enum):
  SPIRAL = auto()
  CONES = auto()

  @classmethod
  def from_index(cls):
    return list(cls)[index]

class Trajectory:
  def __init__(self, file_path) -> None:
      with open(file_path, 'rb') as file:
        version = np.fromfile(file, 'i', 1)[0]
        self.type = TrajectoryType.from_index(np.fromfile(file, 'i', 1)[0])
        num_dimensions = np.fromfile(file, 'i', 1)[0]
        self.image_dimensions = np.fromfile(file, 'i', num_dimensions)
        self.spatial_resolution = np.fromfile(file, 'f', num_dimensions)
        self.field_of_view = np.fromfile(file, 'f', num_dimensions)
        num_readouts = np.fromfile(file, 'i', 1)[0]
        num_bases = np.fromfile(file, 'i', 1)[0]
        max_gradient_amplitude = np.fromfile(file, 'f', 1)[0]
        max_readout_gradient_amplitude = np.fromfile(file, 'f', 1)[0]
        max_slew_rate = np.fromfile(file, 'f', 1)[0]
        num_waveform_points = np.fromfile(file, 'i', 1)[0]
        num_readout_points = np.fromfile(file, 'i', 1)[0]
        sampling_interval = np.fromfile(file, 'f', 1)[0]
        storage = np.fromfile(file, 'i', 1)[0]

        self.gradient_waveforms = read_array(file, 'f').reshape([num_readouts, 3, num_waveform_points])
        gradient_waveform_instructions = read_array(file, 'h')
        self.k_space_coordinates = read_array(file, 'f').reshape([num_readouts, 3, num_readout_points])
        self.density_compensation = read_array(file,  'f').reshape([num_readouts, num_readout_points])


class Cones:
  def __init__(self, file) -> None:
    cones_version = np.fromfile(file, 'i', 1)[0]

    rotatable = np.fromfile(file, 'i', 1)[0]
    interconeCompensation = np.fromfile( file, 'i', 1)[0]
    coneAngles = read_array('f')
    coneAngleDensityCompensation = read_array('f')
    basisConeAngles = read_array('f')
    numBasisReadoutPoints = read_array('i')
    numBasisWaveformPoints = read_array('i')
    basisGradientWaveforms = read_array('f')
    basisKspaceCoordinates = read_array('f')
