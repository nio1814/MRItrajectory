# This source code is under a BSD 3-Clause License.
# See LICENSE for more information.

# To distribute this file, substitute the full license for the above reference.
from enum import Enum
import os
import sys
from PySide2 import QtWidgets, QtUiTools
from PySide2.QtCore import QFile

from cmri.designer.trajectory.generator import QGenerator
from cmri.trajectory.types import TrajectoryType


class Designer(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        
        current_directory = os.path.abspath(os.path.dirname(__file__))
        loader = QtUiTools.QUiLoader()
        file_path_ui = os.path.join(current_directory, 'mainwindow.ui')
        if not os.path.exists(file_path_ui):
            raise FileExistsError(file_path_ui)
        file_ui = QFile(file_path_ui)
        file_ui.open(QFile.ReadOnly)
        if not file_ui:
            raise RuntimeError(f'Failed to open {file_path_ui}')
        self.ui = loader.load(file_ui, parent)
        file_ui.close()

        self.setCentralWidget(self.ui)
        self.generator = QGenerator()
        self.auto_update = True

        self.ui.trajectoryComboBox.currentIndexChanged[int].connect(self.set_trajectory_type)

        for trajectory_type in TrajectoryType:
            self.ui.trajectoryComboBox.addItem(trajectory_type.value, trajectory_type);
        
        self.generator.generate()

    def set_trajectory_type(self, trajectory_type):
        if isinstance(trajectory_type, int):
            trajectory_type = list(TrajectoryType)[trajectory_type]
        self.generator.trajectory_type = trajectory_type


if __name__ == "__main__":
    application = QtWidgets.QApplication()
    designer = Designer()
    designer.show()
    sys.exit(application.exec_())