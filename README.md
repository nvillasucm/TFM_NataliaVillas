# TFM: Calibration and Dosimetry in Nuclear Medicine (Lu-177)

This repository contains the code and resources developed for the Master's Thesis (TFM) by **Natalia Villas**, focused on quantitative SPECT/CT imaging and dosimetry in Lutetium-177 (Lu-177) therapies.
The project provides a graphical tool for phantom-based calibration, enabling accurate estimation of Recovery Coefficients (RC) and activity quantification from clinical imaging data.


## 📌 Project Overview

Accurate dosimetry in nuclear medicine requires reliable quantification of SPECT images. This project addresses that challenge by implementing a digital phantom-based calibration workflow integrated into an interactive GUI.
It includes:
- **Digital Phantom Modeling**:
    - Parametric 2D phantom with hot spheres and background ROIs
    - Adjustable position and rotation for alignment with imaging data
- **SPECT/CT Image Processing**:
    - Loading of DICOM CT series and SPECT volumes
    - Resampling of SPECT to CT spatial grid (XY alignment)
    - Interactive overlay visualization
- **Quantitative Analysis**:
    - Activity estimation in inserts and background regions
    - Automatic detection of the slice with maximum activity
    - Subpixel ROI masking for improved accuracy
- **Recovery Coefficient (RC) Calibration**:
    - Computation of Contrast (Q) and RC
    - Generation of 2D and 3D calibration curves
- **Uncertainity Estimation**:
    - Geometric uncertainty (phantom misalignment sensitivity)
    - Statistical uncertainty (background variability)
    - Combined uncertainty propagation
- **Export of Result**:
    - Graphs (PNG)
    - Numerical results (CSV)

## 🖥️ Graphical User Interface (GUI)
The application is built using Tkinter and designed for clinical/research usability.
Main Capabilities
- Load CT and SPECT DICOM data
- Navigate slices independently (CT vs SPECT)
- Overlay SPECT on CT with transparency
- Enable and manipulate a digital phantom
  - Rotation (°)
  - Translation (pixels)
- Input activity and volume values
- Run:
  - RC/Q calibration

## 📂 Repository Structure


└── Interfaz_calibracion_TFM.ipynb   -->     Main GUI application (Tkinter-based) - Jupyter

└── Interfaz_calibracion_TFM_python.py     -->      Main GUI application (Tkinter-based) - Python File

└── *.dcm       -->          Example SPECT DICOM data

└── README.md

## 🛠 Features

### 1. SPECT/CT Visualization
- Overlay of SPECT functional data on CT anatomical images
- Independent slice navigation for flexible alignment
- Intensity normalization using percentile scaling
- Background masking for improved visualization

### 2. Digital Phantom
- 6 hot spheres:
    - Diameters: 10, 13, 17, 22, 28, 37 mm
- Background ROIs distributed within phantom contour
- Geometry defined in mm and converted to pixels
- Subpixel ROI masks (3×3 sampling) for precision

### 3. Quantitative Analysis
- Mean activity inside:
    - Hot inserts
    - Background regions
- Activity normalized by volume
- Automatic selection of optimal slice (maximum activity)

### 4. Recovery Coefficient (RC) Analysis

### Contrast (Q)
$$
Q = \frac{\left(\frac{A_{insert}}{A_{bg}}\right) - 1}{\left(\frac{a_H}{a_B}\right) - 1} \times 100
$$

### Recovery Coefficient (RC)
$$
RC = Q \cdot \left(1 - \frac{1}{R}\right) + \frac{1}{R}
$$

where:

$$
R = \frac{a_H}{a_B}
$$

### 5. Calibration Curves
- 2D curves (standard)
- 3D curves: $Q^{3/2}$, $RC^{3/2}$
- Error bars include:
    - Systematic (geometric)
    - Statistical uncertainty

## 6. Uncertainty Analysis
Geometric (Systematic)
  - Phantom displacement: ±0.5 pixels
  - Evaluated over a 3×3 grid
Statistical
  - Derived from variability in background ROIs
Total Uncertainty

$$
\sigma_{total} = \sqrt{\sigma_{geom}^2 + \sigma_{stat}^2}
$$

## 🚀 Getting Started

### Prerequisites
To run the interface, you will need Python 3.x and the following libraries:
```bash
pip install numpy matplotlib pydicom SimpleITK opencv-python napari
```

### Running the Interface
1. Open `interfaz.ipynb` in a Jupyter environment (VS Code, JupyterLab, etc.).
2. Run all cells to launch the `tkinter` window.
3. The GUI will launch automatically
4. Select your CT folder and SPECT `.dcm` file to begin the analysis.

### Workflow
1. Load a CT DICOM series
2. Load a SPECT DICOM file
3. Adjust slice positions
4. Enable and align the digital phantom
5. Enter:
    - Hot activity & volume
    - Background activity & volume
6. Run:
    - Q and RC computation
7. Save results (optional)


## 🧪 Methodology
- SPECT is resampled to CT grid (XY only) for spatial consistency
- Z-axis remains independent → flexible slice matching
- ROIs are computed using subpixel masks
- Slice with maximum total activity is automatically selected
- Follows EANM recommendations for quantitative SPECT imaging

## ✍️ Author
**Natalia Villas**  
Master's Thesis (TFM)  
*Keywords: Nuclear Medicine, Lu-177, SPECT/CT, Dosimetry, Medical Imaging, Python, Quantification, Phantom Calibration.* 
