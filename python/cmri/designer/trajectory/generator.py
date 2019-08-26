from enum import Enum

from PySide2.QtCore import QObject

# from cmri.trajectory.generate_spiral

class TrajectoryType(Enum):
    SPIRAL = 'Spiral'
    STACK_OF_SPIRALS = 'Stack of Spirals'


class Generator:
    def __init__(self):
        self.field_of_view = [28] * 3
        self.spatial_resolutoin = [2] * 3
        self.trajectory_type = TrajectoryType.SPIRAL

    def generate(self):
        if self.trajectory_type == TrajectoryType.SPIRAL:
            pass


class QGenerator(Generator, QObject):
    def __init__(self, parent=None):
        super(QObject, self).__init__()
