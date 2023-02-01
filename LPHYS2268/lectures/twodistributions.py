# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 11:42:28 2020

@author: massonnetf

@purpose: illustration of marginal distribution
"""

# Imports
# -------
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy import stats

# Numerical methods for solving the system 
# of ordinary differential equations (ODEs)
from numerical_methods import *

# Seed for reproducibility
np.random.seed(0)

M = 10000 # State dimension

fig = plt.figure(figsize = (4, 2), dpi = 300)

x1 = np.exp(0.5 * np.random.randn(M))

x2 = np.exp(np.random.randn(M))

bins = np.arange(-3.0, 10.0, 0.1)
plt.hist(x1, bins, alpha = 0.7, density = True, 
         label = " Forecast distribution")
plt.hist(x2, bins, alpha = 0.7, density = True, 
         label = " Climatological distribution")
plt.xlim(0.0, 6)
plt.xlabel("Precipitation in Louvain-la-Neuve tomorrow [mm]")
plt.ylabel("pdf")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("./fig6.png")
