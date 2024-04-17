# Soft Photon Heating

This Jupyter notebook allows one to study the effects of soft photon heating as described in 2024.xxxx.

## Installation
Simply download the three files, place them in the directory of your choice, and run the Soft-Photon-Heating.ipynb. 

## File contents

### Soft-Photon-Heating.ipynb
Contains the main elements necessary to investigate the response of the matter temperature (T_M), ionization fraction (X_e), and spectral distortion (\Delta n) to the injection of soft photon backgrounds.

To use: modify relevant sections in the 'User Input Section' cell, and then run rest of the cells.
Plots of T_M and X_e compare against pre-tabulated solutions from CosmoRec (vanilla LCDM).
The cosmos.py contains details of most of the function calls used by the solver.

### cosmos.py
Contains various fundamental parameters and constants, as well as commonly used functions which are called extensively by Soft-Photon-Heating.ipynb. The code operates using natural units, and the natural unit associated with each object in this file are labelled using square [] brackets. 

### CosmoRec.Sol.dat
Contains highly accurate data on the matter temperature and ionization fraction in the case of no soft photon injections, at redshifts 100 < z < 8500. The data is computed using CosmoRec ([Chluba and Thomas, 2011](https://inspirehep.net/literature/873187)), and is used in the Jupyter notebook for comparing the soft photon solutions against the vanilla LCDM evolution of T_M and X_e.
