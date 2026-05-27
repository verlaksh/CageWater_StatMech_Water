import os ; import sys ; import matplotlib.pyplot as plt ; from math import *; import numpy as np
from scipy.special import erf ; from scipy.special import expi ; from scipy.integrate import quad


#Default parameters

# r_HB and eps_HB real units
eps_HB_real= 19.200             # HB potential in kJ/mol
r_HB_real  = 0.277               # HB distance in nm
T_real     = 298.15                # T in K
P_real     = 101325             # P in Pascal

# Constants
NA         = 6.022*10**23       # Avogadro number in per mole
kB         = 1.381*10**-23      # Boltzmann constant in J/K
R          = NA*kB              # Gas Constant J/molK
nm_to_m    = 10**-9             # nm to meter conversion

## parameters Parameters in reduced units
a          =  0.33300            # van der Waal's constant
eps_HB     =  1.00000            # HB potential
eps_LJ     =  0.53690            # LJ potential
eps_c      = -0.07000            # Cooperative constant for cage water
KT_HB      =  6.00000            # Spring constant for cage and HB 
KT_LJ      =  6.00000            # Spring constant for LJ 
Kx         =  5.25490            # Tetrahedral spring constant
Kth        = 46.72980            # Radial spring constant
r_HB       =  1.00000            # HB distance
sigma_LJ   =  1.00000            # Average size in Ih ice
r_d        =  0.05000            # hardcore size of water # Size of water molecule
r_mS       =  0.07200            # Max bond length cage
r_mHB      =  0.18000            # Max bond length HB
r_mLJ      =  0.27290            # Max bond length LJ
xHB        =  0.72350            # Perturbation for HB
xLJ        =  1.28510            # Perturbation for LJ


#for interactive runs
#19.2    0.333   0.5369  -0.07   6       6       5.2549  46.7298 0.05    0.072   0.18    0.2729  0.7235  1.2851
