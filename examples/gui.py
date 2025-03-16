from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QPushButton, QFileDialog, QWidget
from PyQt5.QtWebEngineWidgets import QWebEngineView
import os

class MolecularViewerApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molecular Viewer")
        self.setGeometry(100, 100, 800, 600)

        # Main widget
        widget = QWidget()
        layout = QVBoxLayout()
        widget.setLayout(layout)

        # 3Dmol.js Viewer
        self.viewer = QWebEngineView()
        self.viewer.setHtml(self._generate_html())
        layout.addWidget(self.viewer)

        # Load File Button
        self.load_button = QPushButton("Load Molecule")
        self.load_button.clicked.connect(self.load_molecule)
        layout.addWidget(self.load_button)

        self.setCentralWidget(widget)

    def _generate_html(self, molecule_data=""):
        return f"""
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://3dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        </head>
        <body>
            <div style="width: 100%; height: 100%;" id="viewer"></div>
            <script>
                let viewer = $3Dmol.createViewer("viewer", {{defaultcolors: $3Dmol.rasmolElementColors}});
                if ('{molecule_data}') {{
                    viewer.addModel(`{molecule_data}`, "xyz");
                    viewer.setStyle({{}}, {{stick: {{}}, sphere: {{}}}});
                    viewer.zoomTo();
                    viewer.render();
                }}
            </script>
        </body>
        </html>
        """

    def load_molecule(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Open Molecule File", "", "XYZ Files (*.xyz);;PDB Files (*.pdb)")
        if file_path:
            with open(file_path, "r") as file:
                molecule_data = file.read()
                self.viewer.setHtml(self._generate_html(molecule_data))

if __name__ == "__main__":
    app = QApplication([])
    window = MolecularViewerApp()
    window.show()
    app.exec_()
