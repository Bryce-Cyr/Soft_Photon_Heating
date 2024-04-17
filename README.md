# Soft Photon Heating
This Jupyter notebook allows one to study the effects of soft photon heating as described in 2024.xxxx.

## Installation
Simply download the three files, place them in the directory of your choice, and run the Soft-Photon-Heating.ipynb. 

## File contents

### Soft-Photon-Heating.ipynb
Contains the main elements necessary to investigate the response of the matter temperature ($T_M$), ionization fraction ($X_e$), and spectral distortion ($\Delta n$) to the injection of soft photon backgrounds. Regions marked "Switches" can be set in any way the user desires. Regions marked "Select one" means to just turn one of the options on at a time. 

To use: modify relevant sections in the 'User Input Section' cell, and then run rest of the cells.

Plotting code can be used as is, or by modifying switches present in those cells. T_M and X_e compare against pre-tabulated (vanilla LCDM) solutions from CosmoRec.

The cosmos.py contains details of most of the function calls used by the solver.

### cosmos.py
Contains various fundamental parameters and constants, as well as commonly used functions which are called extensively by Soft-Photon-Heating.ipynb. The code operates using natural units, and the natural unit associated with each object in this file are labelled using square [] brackets. Modification of this file is not usually necessary.

### CosmoRec.Sol.dat
Contains highly accurate data on the matter temperature and ionization fraction in the case of no soft photon injections, at redshifts 100 < z < 8500. The data is computed using CosmoRec ([Chluba and Thomas, 2011](https://inspirehep.net/literature/873187)), and is used in the Jupyter notebook for comparing the soft photon solutions against the vanilla LCDM evolution of T_M and X_e.

## Warnings and limitations
- Most importantly, this code does not contain a reionization module. That means the solutions computed are fully accurate over 30 < z < 1500. Please interface this code with your favourite reionization module to properly evolve to low redshifts.
- The runtime
