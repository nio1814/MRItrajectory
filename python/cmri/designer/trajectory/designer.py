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
        self.link_xy = False
        
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

        self.ui.trajectoryComboBox.currentIndexChanged[int].connect(self.set_trajectory_type)

        for trajectory_type in TrajectoryType:
            self.ui.trajectoryComboBox.addItem(trajectory_type.value, trajectory_type);
        
        self.field_of_view_spin_boxes = [
            self.ui.fieldOfViewXSpinBox, 
            self.ui.fieldOfViewYSpinBox, 
            self.ui.fieldOfViewZSpinBox
        ]
        for axis, spin_box in enumerate(self.field_of_view_spin_boxes):
            spin_box.valueChanged.connect(lambda fov, dimension=axis: self.set_field_of_view(fov, dimension))

        self.field_of_view_sliders = [
            self.ui.fieldOfViewXSlider, 
            self.ui.fieldOfViewYSlider, 
            self.ui.fieldOfViewZSlider
        ]
        for axis, slider in enumerate(self.field_of_view_sliders):
            slider.valueChanged.connect(lambda fov, dimension=axis: self.set_field_of_view(fov, dimension))

        self.ui.generatePushButton.clicked.connect(self.generator.generate)

        field_of_view_min = 140
        field_of_view_max = 400
        spatial_resolution_min = .5
        spatial_resolution_max = 4
        for axis in range(3):
            self.field_of_view_sliders[axis].setMinimum(field_of_view_min)
            self.field_of_view_sliders[axis].setMaximum(field_of_view_max)
            self.field_of_view_spin_boxes[axis].setMinimum(field_of_view_min)
            self.field_of_view_spin_boxes[axis].setMaximum(field_of_view_max)
            self.set_field_of_view(280, axis)

        self.set_trajectory_type(TrajectoryType.SPIRAL)

    def set_field_of_view(self, field_of_view, axis):
        self.generator.field_of_view[axis] = field_of_view

        spin_box = self.field_of_view_spin_boxes[axis]
        slider = self.field_of_view_sliders[axis]
        for element in [spin_box, slider]:
            if element.value != field_of_view:
                element.setValue(field_of_view)
        
        if self.link_xy:
            other_axis = 1 - axis
            other_elments = [
                self.field_of_view_spin_boxes[other_axis],
                self.field_of_view_spin_boxes[other_axis]
            ]
            if any([element.value() != field_of_view for element in other_elments]):
                self.set_field_of_view(field_of_view, other_axis)

    def set_trajectory_type(self, trajectory_type):
        if isinstance(trajectory_type, int):
            trajectory_type = list(TrajectoryType)[trajectory_type]
        self.generator.trajectory_type = trajectory_type

        if trajectory_type == TrajectoryType.SPIRAL:
            self.link_xy = True


if __name__ == "__main__":
    application = QtWidgets.QApplication()
    designer = Designer()
    designer.show()
    sys.exit(application.exec_())
