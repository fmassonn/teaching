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

[ax.plot3D(t * time_unit, X[jMC, jM, :], np.zeros(len(t)), lw = 0.01, \
     color = "grey") \
     for jMC in range(nMC)]
# Plot truth
ax.plot3D(t * time_unit, X[jtruth, jM, :], np.zeros(len(t)), lw = 1, \
          color = [0.0, 0.5, 0.0])

ax.plot3D((0.0, nd_show), (0.0, 0.0), (0.0, 0.0), "k--")
# Estimating density for several times
x_pdf = np.linspace(-2.0, 4.0, 1000)
for day in np.arange(0, nd_show, 20):
    kernel = stats.gaussian_kde(X[:, 0, int(day / (time_unit * dt))])
    pdf    = kernel(x_pdf).T
    ax.plot3D(np.full(len(pdf), day), x_pdf, pdf, "orange", lw = 1)

plt.grid()
plt.xlim(0, t[-1] * time_unit)
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





