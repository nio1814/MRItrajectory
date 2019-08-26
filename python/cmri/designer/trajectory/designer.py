import os
import sys
from PySide2 import QtWidgets, QtUiTools
from PySide2.QtCore import QFile

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
        

if __name__ == "__main__":
    application = QtWidgets.QApplication()
    designer = Designer()
    designer.show()
    sys.exit(application.exec_())