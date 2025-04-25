# PTSFit - Pressure, Temperature from Sequential Fits
  
## 1. Overview:

This program is designed to perform Le-Bail fits on a large number of sequential X-ray diffractograms in order to determine unit cell parameters and volume.
A .cif file may be used as a starting model, or a user can input unit cell and symmetry parameters.
The primary use case is the determination of unit cell volume of a pressure/temperature calibrant present from in-situ X-ray diffraction data.
This provides a quick and efficient method for determination of pressure and temperature conditions from the equation of state of the calibrant phase.
The program uses a BM EoS function to determine P and T from the refined crystallographic volume of the calibrant phase. Built-in EoS parameters are provided for common calibrants or may be input by the user.
  
Rudimentary plotting functions are provided in the GUI, these are meant as a diagnostic tool rather than for producing publication quality figures.
The results of the sequential refinement, PVT parameters, and derived thermodynamic parameters associated with the BM equation of state can be saved to a .csv file.
  
Analysis of more complex, lower symmetry crystal phases is possible but the peak overlap common in lower symmetry systems is likely to cause instabilities in Le-Bail fitting techniqes. For these cases dedicated Rietveld refinement software such as FullProf or GSAS / GSAS II is recommended.  

## 2. Citing:
TBA

## 3. Installation:

Extract the released .zip archive and run .exe file

## 4. User defined squation of state parameters:
  
Currently, the PVT.py file contains a single class for determination of P,T from a BM EoS. Other EoS functions can be added as seperate classes. The EoS_dictionaries file contains dictionaries of in-built EoS parameters for commonly used calibrant phases. Users can add parameters here for other phases.
If you would like to add an in-built EoS dictionary to the .exe version the relavent files are localted in _internal.

## 5. References:
PVT calculation based on Sokolova et al: https://doi.org/10.1016/j.cageo.2016.06.002  
Temperature calculation derived from McHardy et al: https://doi.org/10.1080/08957959.2023.2187294  
Indexing via xrayutilities from Kreigner et al. https://doi.org/10.1107/S0021889813017214  
