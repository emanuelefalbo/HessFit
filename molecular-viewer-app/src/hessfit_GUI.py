import tkinter as tk
from tkinter import filedialog, messagebox
import os
import json
import py3Dmol

class MolecularViewerApp:
    def __init__(self, master):
        self.master = master
        master.title("Molecular Viewer")

        self.viewer = py3Dmol.view(width=800, height=400)
        self.viewer.show()

        self.load_button = tk.Button(master, text="Load Molecule", command=self.load_molecule)
        self.load_button.pack()

    def load_molecule(self):
        json_file = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
        if json_file:
            self.load_molecule_from_json(json_file)

    def load_molecule_from_json(self, json_file):
        if not os.path.exists(json_file):
            messagebox.showerror("Error", "JSON file not found!")
            return

        try:
            with open(json_file, "r") as f:
                data = json.load(f)

            molecule_file = data.get("molecule")
            if not molecule_file or not os.path.exists(molecule_file):
                messagebox.showerror("Error", "Molecule file not found in JSON or invalid path.")
                return

            with open(molecule_file, "r") as mol_file:
                molecule_data = mol_file.read()

            file_extension = os.path.splitext(molecule_file)[-1].lower()
            self.viewer.removeAllModels()
            if file_extension == ".xyz":
                self.viewer.addModel(molecule_data, "xyz")
            elif file_extension == ".pdb":
                self.viewer.addModel(molecule_data, "pdb")
            else:
                messagebox.showerror("Error", "Unsupported file format.")
                return

            self.viewer.setStyle({}, {"stick": {}, "sphere": {}})
            self.viewer.zoomTo()
            self.viewer.show()
        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    root = tk.Tk()
    app = MolecularViewerApp(root)
    root.mainloop()