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