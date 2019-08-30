from PySide2.QtCore import QObject, Signal

from cmri.trajectory.generator import Generator


class QGenerator(Generator, QObject):
    generated = Signal(dict)

    def __init__(self):
        super().__init__()
        QObject.__init__(self)
        self.auto_update = True

    def generate(self):
        trajectory = super().generate()
        if trajectory:
            self.generated.emit(trajectory)

    def set_axis_field_of_view(self, field_of_view, axis):
        field_of_view_new = list(self.field_of_view)
        field_of_view_new[axis] = field_of_view

        num_dimensions = self.trajectory_type.num_dimensions()
        if self.field_of_view[:num_dimensions] != field_of_view_new[:num_dimensions]:
            self.field_of_view = field_of_view_new
            if self.auto_update:
                self.generate()

    def set_readout_duration(self, readout_duration):
        if self.readout_duration != readout_duration:
            self.readout_duration = readout_duration
            if self.auto_update:
                self.generate()