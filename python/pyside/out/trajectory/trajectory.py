from enum import Enum

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
  SPIRAL = 'spiral'
  RADIAL = 'radial'
  RADIAL3D = 'radial3d'
  CONES = 'cones'

  @classmethod
  def from_index(cls, index):
    return list(cls)[index]

class Trajectory:
  def __init__(self, file_path):
      with open(file_path, 'rb') as file:
        version = np.fromfile(file, 'i', 1)[0]
        self.type = TrajectoryType.from_index(np.fromfile(file, 'i', 1)[0])
        num_dimensions = np.fromfile(file, 'i', 1)[0]
        self.image_dimensions = np.fromfile(file, 'i', num_dimensions)
        self.spatial_resolution = np.fromfile(file, 'f', num_dimensions)
        self.field_of_view = np.fromfile(file, 'f', num_dimensions)
        num_readouts = np.fromfile(file, 'i', 1)[0]
        num_bases = np.fromfile(file, 'i', 1)[0]
        self.max_gradient_amplitude = np.fromfile(file, 'f', 1)[0]
        self.max_readout_gradient_amplitude = np.fromfile(file, 'f', 1)[0]
        self.max_slew_rate = np.fromfile(file, 'f', 1)[0]
        num_waveform_points = np.fromfile(file, 'i', 1)[0]
        num_readout_points = np.fromfile(file, 'i', 1)[0]
        self.sampling_interval = np.fromfile(file, 'f', 1)[0]
        storage = np.fromfile(file, 'i', 1)[0]

        self.gradient_waveforms = read_array(file, 'f').reshape([num_readouts, 3, num_waveform_points])
        gradient_waveform_instructions = read_array(file, 'h')
        self.k_space_coordinates = read_array(file, 'f').reshape([num_readouts, 3, num_readout_points])
        self.density_compensation = read_array(file,  'f').reshape([num_readouts, num_readout_points])

        if self.type == TrajectoryType.CONES:
          self.cones = Cones(file)


class Cones:
  def __init__(self, file) -> None:
    version = np.fromfile(file, 'i', 1)[0]

    rotatable = np.fromfile(file, 'i', 1)[0]
    self.intercone_compensation = np.fromfile( file, 'i', 1)[0]
    self.cone_angles = read_array(file, 'f')
    self.cone_angle_density_compensation = read_array(file, 'f')
    self.basis_cone_angles = read_array(file, 'f')
    self.num_basis_readout_points = read_array(file, 'i')
    self.num_basis_waveform_points = read_array(file, 'i')
    self.basis_gradient_waveforms = read_array(file, 'f')
    self.basis_k_space_coordinates = read_array(file, 'f')
