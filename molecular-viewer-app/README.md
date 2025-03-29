# Molecular Viewer App

This project is a molecular viewer application that allows users to load and visualize molecular data using two different GUI frameworks: Tkinter and PyQt. The application supports loading molecular files in formats such as XYZ and PDB, and it utilizes 3Dmol.js for rendering molecular structures.

## Project Structure

```
molecular-viewer-app
├── src
│   ├── hessfit_GUI.py       # Main logic for loading and rendering molecular data using Tkinter
│   ├── gui.py               # GUI setup using PyQt with a 3Dmol.js viewer
│   └── utils
│       └── __init__.py      # Package initializer for utility functions
├── requirements.txt          # List of dependencies required for the project
└── README.md                 # Documentation for the project
```

## Setup Instructions

1. **Clone the repository**:
   ```
   git clone <repository-url>
   cd molecular-viewer-app
   ```

2. **Install dependencies**:
   Make sure you have Python installed, then run:
   ```
   pip install -r requirements.txt
   ```

## Usage

- To run the Tkinter version of the molecular viewer, execute:
  ```
  python src/hessfit_GUI.py
  ```

- To run the PyQt version of the molecular viewer, execute:
  ```
  python src/gui.py
  ```

## Dependencies

This project requires the following Python packages:
- Tkinter
- PyQt5
- py3Dmol
- Any other necessary packages for molecular visualization

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue for any enhancements or bug fixes.

## License

This project is licensed under the MIT License. See the LICENSE file for more details.