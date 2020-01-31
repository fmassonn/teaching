# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:10:05 2020

@author: massonnetf

@reference: Lorenz, Edward N. “Irregularity: A Fundamental Property of the 
            Atmosphere.” Tellus A 36A, no. 2 (1984): 98–110. 
            https://doi.org/10.1111/j.1600-0870.1984.tb00230.x.

"""

# Imports
# -------
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy import stats

# Numerical methods for solving the system 
# of ordinary differential equations (ODEs)
from numerical_methods import *

# Seed for reproducibility
np.random.seed(0)

# Model parameters
# State vector dimension
M = 3
# Index of observable
jM = np.random.randint(M)

a = 0.25
b = 4.0
F = 8.0
G = 1.0

# Run parameters
t0 = 0.0 # Initial time
# One time unit is 5 days
time_unit = 5.0
# Time step expressed in time units
dt = 1 / 30
# Integration time (in days), then conversion from days to time steps
nd = 500
# Nb days to show
nd_show = 90 
N = int(nd / (time_unit * dt))
# The  time axis
t = np.array([t0 + n * dt for n in range(N)])

# Number of Monte-Carlo integrations
nMC = 1000
std = 1e-1 # standard deviation of Gaussian in perturbation
# Index of truth
jtruth = np.random.randint(nMC) 

# The initial condition
X0 = np.full(M, 1) 

# The system to be solved
#
# dX/dt = - Y^2 - Z^2 - aX + aF
# dY/dt =   XY  - bXZ - Y  + G
# dZ/dt =  bXY  + XZ  - Z

# Right-hand side of ODE as a function
# X is a column np array of size M
def f(t, M, X):
        
    value = np.full(M, np.nan)
    
    value[0] = - X[1] ** 2 - X[2] ** 2 - a * X[0] + a * F
    value[1] = X[0] * X[1] - b * X[0] * X[2] - X[1] + G
    value[2] = b * X[0] * X[1] + X[0] * X[2] - X[2]
    
    return value

X = np.full((nMC, M, N), np.nan)



# Start Monte-Carlo integrations
for jMC in np.arange(nMC):
    print(jMC)
    # Start the run
    # -------------
    X[jMC, :, 0] = X0 + std * np.random.randn(M)
    
    for jN in np.arange(1, N):
        # Runge-Kutta 4th order method
        X[jMC, :, jN] = rk4vec(t[jN - 1], M, X[jMC, :, jN - 1], dt, f)


# Plot
# ----
fig, ax = plt.subplots(figsize = (6, 5), dpi = 300, \
                       constrained_layout = True)
ax = plt.axes(projection = '3d')
ax.view_init(50, -45)
jd_show = int(nd_show / (time_unit * dt))

has_legend = False
for jMC in range(nMC):
    if not has_legend:
        label = "Realization"
        has_legend = True
    else:
        label = None
      
    ax.plot3D(t[:jd_show] * time_unit, X[jMC, jM, :jd_show], \
           np.zeros(len(t[:jd_show])), lw = 0.01, \
    color = "grey", label = label)

# Plot truth
ax.plot3D(t[:jd_show] * time_unit, X[jtruth, jM, :jd_show], \
          np.zeros(len(t[:jd_show])), lw = 1, \
          color = [0.0, 0.5, 0.0], label = "Truth")

ax.plot3D((0.0, nd_show), (0.0, 0.0), (0.0, 0.0), "k--")
# Estimating density for several times
x_pdf = np.linspace(-2.0, 4.0, 1000)
has_legend = False
for day in np.arange(0, nd_show, 30):
    kernel = stats.gaussian_kde(X[:, jM, int(day / (time_unit * dt))])
    pdf    = kernel(x_pdf).T
    if not has_legend:
        label = "Ensemble PDF"
        has_legend = True
    else:
        label = None
        
    ax.plot3D(np.full(len(pdf), day), x_pdf, pdf, "orange", lw = 1, \
              label = label)

# Fit climatology
kernel_clim = stats.gaussian_kde(X[jtruth, jM, :])
pdf_clim = kernel_clim(x_pdf).T
ax.plot3D(np.full(len(pdf_clim), 0), x_pdf, pdf_clim, \
          color = [1.0, 0.5, 0.5], lw = 2, label = "Climatology PDF")
plt.legend()
plt.grid()
plt.xlim(0, t[jd_show] * time_unit)
plt.xlabel("Days")
plt.ylabel("Observable")
plt.ylim(-2.0, 4.0)   
ax.set_zlim(0.1, 4.0)
ax.set_zlabel("PDF")
plt.tight_layout(pad = 5.0)
plt.savefig("./fig001.png")

## 3D
#fig, ax = plt.subplots(figsize = (6, 5), dpi = 100, constrained_layout=True)
#ax = plt.axes(projection='3d')
#ax.view_init(50, -45)
#[ax.plot3D(X[jMC, 0, :200], X[jMC, 1, :200], X[jMC, 2, :200], "k.", ms = 0.5, \
#           color = [0.5, 0.5, 0.5]) for jMC in range(nMC)]
#ax.plot3D(X[jtruth, 0, :], X[jtruth, 1, :], X[jtruth, 2, :], ls = "-", lw = 1, \
#           color = [0.0, 0.5, 0.0])
#plt.tight_layout(pad = 5.0)
#plt.savefig("./fig002.png")





