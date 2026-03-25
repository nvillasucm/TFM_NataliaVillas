# TFM: Calibration and Dosimetry in Nuclear Medicine (Lu-177)

This repository contains the code and resources developed for the Master's Thesis (TFM) by **Natalia Villas**. The project focuses on the calibration and calculation of Recovery Coefficients (RC) for Lu-177 treatments using SPECT/CT imaging and digital phantom simulations.

## 📌 Project Overview

The main objective of this work is to provide a tool for the analysis of nuclear medicine images (SPECT/CT) to improve dosimetry accuracy in Lutetium-177 therapies. It includes:
- **Digital Phantom Generation**: Modeling of standardized phantoms with hot inserts and background spheres.
- **Image Processing**: Resampling and overlaying SPECT data onto CT grids for precise localization.
- **Recovery Coefficient (RC) Calculation**: Tools to calculate activity concentration and RC in 2D and 3D with statistical and systematic uncertainty estimation.
- **Graphical User Interface (GUI)**: A specialized interface for clinicians and researchers to interact with DICOM series.

## 📂 Repository Structure

- `interfaz.ipynb`: The main application. A Python-based GUI (using `tkinter`) for loading DICOM series, visualizing SPECT/CT overlays, and performing RC calculations.
- `mascaras.ipynb`: Notebook focused on the creation of masks and digital phantoms using `napari` and `OpenCV`.
- `Tomo4FOV_Lu177peak_IRACSC001_DS.dcm`: Sample DICOM file containing Lu-177 SPECT data.

## 🛠 Features

### 1. SPECT/CT Interface
The interface allows users to:
- Load CT series and SPECT files.
- Visualize synchronized slices with adjustable transparency (alpha-blending).
- Calculate activity in specific Regions of Interest (ROIs) corresponding to phantom inserts.
- Estimate uncertainties (systematic and statistical) for the calculated coefficients.

### 2. Recovery Coefficient (RC) Analysis
- **2D/3D Curves**: Generation of RC curves based on sphere diameters.
- **Uncertainty Propagation**: Calculation of errors derived from voxel positioning and background statistics.
- **Exporting Results**: Save graphs as PNG and data tables as CSV.

## 🚀 Getting Started

### Prerequisites
To run the notebooks and the interface, you will need Python 3.x and the following libraries:
```bash
pip install numpy matplotlib pydicom SimpleITK opencv-python napari
```

### Running the Interface
1. Open `interfaz.ipynb` in a Jupyter environment (VS Code, JupyterLab, etc.).
2. Run all cells to launch the `tkinter` window.
3. Select your CT folder and SPECT `.dcm` file to begin the analysis.

## 🧪 Methodology
The project utilizes a 2D/3D digital phantom consisting of:
- **6 Hot Inserts**: Diameters of 10, 13, 17, 22, 28, and 37 mm.
- **Background ROIs**: Used to calculate the target-to-background ratio and statistical noise.

The calculations follow the EANM recommendations for dosimetry in internal radiotherapy.

## ✍️ Author
**Natalia Villas**  
Master's Thesis (TFM)  
*Keywords: Nuclear Medicine, Lu-177, SPECT/CT, Dosimetry, Python, Medical Imaging.* 
