from PySide2.QtCore import QObject

from cmri.trajectory.generator import Generator


class QGenerator(Generator, QObject):
    def __init__(self, parent=None):
        super(QObject, self).__init__()
