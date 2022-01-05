# This source code is under a BSD 3-Clause License.
# See LICENSE for more information.

# To distribute this file, substitute the full license for the above reference.
import os
import sys
from logging import warning

from PySide2 import QtWidgets, QtUiTools
from PySide2.QtCore import QFile
import numpy as np
import pyqtgraph as pg

from cmri.cmri import TrajectoryGenerator, TrajectoryType

TRAJECTORY_NAMES = {TrajectoryType.CONES: 'Cones',
                    TrajectoryType.SPIRAL: 'Spiral'}

def invert_dictionary(data):
    return {value: key for key, value in data.items()}

def name_to_trajectory_type(name):
    return invert_dictionary(TRAJECTORY_NAMES)[name]

class Designer(QtWidgets.QMainWindow):
    READOUT_DURATION_SLIDER_SCALE = 1e4
    READOUT_DURATION_SPIN_BOX_SCALE = 1e3

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
        self.generator = TrajectoryGenerator()

        self.gradients_plot = pg.PlotWidget(name='Gradients')
        self.gradients_plot.setXRange(0, self.ui.readoutDurationSpinBox.maximum())
        self.gradients_plot.setYRange(-4, 4)
        self.ui.plotGridLayout.addWidget(self.gradients_plot, 0, 1);
        self.ui.trajectoryComboBox.currentIndexChanged[str].connect(self.set_trajectory_type)

        for trajectory_type, name in TRAJECTORY_NAMES.items():
            self.ui.trajectoryComboBox.addItem(name, trajectory_type)

        self.field_of_view_spin_boxes = [
            self.ui.fieldOfViewXSpinBox,
            self.ui.fieldOfViewYSpinBox,
            self.ui.fieldOfViewZSpinBox
        ]
        for axis, spin_box in enumerate(self.field_of_view_spin_boxes):
            spin_box.editingFinished.connect(lambda dimension=axis: self.set_field_of_view(None, dimension))

        self.field_of_view_sliders = [
            self.ui.fieldOfViewXSlider,
            self.ui.fieldOfViewYSlider,
            self.ui.fieldOfViewZSlider
        ]
        for axis, slider in enumerate(self.field_of_view_sliders):
            slider.valueChanged.connect(lambda fov, dimension=axis: self.set_field_of_view(fov, dimension))

        self.ui.readoutDurationSlider.sliderMoved.connect(lambda index: self.set_readout_duration(index / self.READOUT_DURATION_SLIDER_SCALE))
        self.ui.readoutDurationSpinBox.editingFinished.connect(lambda: self.set_readout_duration(None))

        self.ui.generatePushButton.clicked.connect(self._generate)

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
        self.set_readout_duration(4e-3)
        self.generator.generate()

    def set_field_of_view(self, field_of_view, axis, source=None):
        if axis >= self.generator.numDimensions():
            warning(f'Attempting to set {field_of_view} fov on axis {axis} for {self.generator.trajectoryType()} '
                    'trajectory')
            return

        if source == 'spin':
            field_of_view = self.field_of_view_spin_boxes[axis].value()

        axes = [axis]
        if self.link_xy and axis in [0, 1]:
            axes.append(1 - axis)
        for axis in axes:
            spin_box = self.field_of_view_spin_boxes[axis]
            slider = self.field_of_view_sliders[axis]
            for element in [spin_box, slider]:
                element.blockSignals(True)
                if element.value() != field_of_view:
                    element.setValue(field_of_view)
                element.blockSignals(False)
            self.generator.setAxisFieldOfView(field_of_view / 10, axis)

    def set_readout_duration(self, readout_duration):
        if readout_duration is None:
            readout_duration = self.ui.readoutDurationSpinBox.value() / self.READOUT_DURATION_SPIN_BOX_SCALE
        self.ui.readoutDurationSlider.setValue(readout_duration * self.READOUT_DURATION_SLIDER_SCALE)
        self.ui.readoutDurationSpinBox.setValue(readout_duration * self.READOUT_DURATION_SPIN_BOX_SCALE)
        self.generator.setReadoutDuration(readout_duration)

    def set_trajectory_type(self, trajectory_type):
        if isinstance(trajectory_type, str):
            trajectory_type = name_to_trajectory_type(trajectory_type)
        self.generator.trajectory_type = trajectory_type

        if trajectory_type == TrajectoryType.SPIRAL:
            self.link_xy = True

    def _generate(self):
        self.generator.generate()
        self.update_plots(self.generator.trajectory())

    def update_plots(self, trajectory):
        readout_index = self.ui.readoutSlider.value()
        self.gradients_plot.clear()
        data = trajectory.gradient_waveforms()[readout_index]
        time = 1e3 * trajectory.samplingInterval * np.arange(data.shape[-1])
        for axis, gradient in enumerate(data):
            curve = self.gradients_plot.plot(x=time, y=gradient)
            color = [255] * 3
            color[axis] = 0
            curve.setPen(color)


if __name__ == "__main__":
    application = QtWidgets.QApplication()
    designer = Designer()
    designer.show()
    sys.exit(application.exec_())
