# PTSeqFit - Pressure, Temperature from Sequential Fits
  
## 1. Overview


This program is designed to perform Le-Bail fits on a large number of sequential X-ray diffractograms in order to determine unit cell parameters and volume.
The primary use case is the determination of unit cell volume of a pressure/temperature calibrant.
This provides a quick and efficient method for determination of pressure and temperature conditions from the equation of state of the calibrant phase.
The program uses a BM EoS function to determine P and T from the refined crystallographic volume of the calibrant phase. Built-in EoS parameters are provided for common calibrants or may be input by the user, or read from appropriate JCPDS files (not yet implemented).  
  
Analysis of more complex, lower symmetry crystal phases is possible but the peak overlap common in lower symmetry systems is likely to cause instabilities in Le-Bail fitting techniqes. For these cases dedicated Rietveld refinement software such as FullProf or GSAS / GSAS II is recommended.  

## 2. Installation

Probably will generate a .exe at some point.

## 3. User defined squation of state
  
Currently, the PVT.py file contains a single class for determination of P,T from a BM EoS. Other EoS functions can be added as seperate classes. The EoS_dictionaries file contains dictionaries of in-built EoS parameters for commonly used calibrant phases. Users can add parameters here for other phases.

