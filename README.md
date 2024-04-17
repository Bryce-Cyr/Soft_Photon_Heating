# Soft Photon Heating
This Jupyter notebook allows one to study the effects of soft photon heating as described in 2024.xxxx.

## Installation
Simply download the three files, place them in the directory of your choice, and run the Soft-Photon-Heating.ipynb. 

## File contents

### Soft-Photon-Heating.ipynb
Contains the main elements necessary to investigate the response of the matter temperature ($T_M$), ionization fraction ($X_e$), and spectral distortion ($\Delta n$) to the injection of soft photon backgrounds. Please see the below for warnings and limitations of the code.

To use: modify relevant sections in the 'User Input Section' cell, and then run rest of the cells. Regions marked "Switches" can be set in any way the user desires. Regions marked "Select one" means to just turn one of the options on at a time. The pre-loaded settings use the source term from Eq. 21 of the paper with $\gamma = 3.6$, $\Delta \rho/\rho = 10^{-6}$, $z_{\rm inj} = 500$, and $E_{\rm cut} = 0.235$ eV. This takes roughly 20 seconds to run.

The plotting code can be used as is, or by modifying switches present in those cells. T_M and X_e compare against pre-tabulated (vanilla LCDM) solutions from CosmoRec.

The cosmos.py contains details of most of the function calls used by the solver.

### cosmos.py
Contains various fundamental parameters and constants, as well as commonly used functions which are called extensively by Soft-Photon-Heating.ipynb. The code operates using natural units, and the natural unit associated with each object in this file are labelled using square [] brackets. Modification of this file is not usually necessary.

### CosmoRec.Sol.dat
Contains highly accurate data on the matter temperature and ionization fraction in the case of no soft photon injections, at redshifts $100 < z < 8500$. The data is computed using CosmoRec ([Chluba and Thomas, 2011](https://inspirehep.net/literature/873187)), and is used in the Jupyter notebook for comparing the soft photon solutions against the vanilla LCDM evolution of T_M and X_e.

## Warnings and limitations
- Most importantly, this code does not contain a reionization module. That means the solutions computed are fully accurate over $30 \lesssim z < 1500$. Please interface this code with your favourite reionization module to properly evolve to low redshifts.
- When running the solver for the first time, it is possible to get overflow warnings from the recombination calculation. This does not affect the results of the code and can be ignored.
- The runtime is relatively insensitive to the photon frequency resolution (set by x_fidel), but quite sensitive to the redshift step sizes (z_step). The default values x_fidel = 5000 and z_step = 0.2 should be able to handle any set of parameters with reasonable energy injections.
- If injecting large amounts of energy ($\delta \rho/\rho \gtrsim 10^{-5}$) around recombination ($900 \lesssim z \lesssim 1200$), the solver might crash. If it does, simply reduce z_step and try again.
