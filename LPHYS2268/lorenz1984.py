# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:10:05 2020

@author: massonnetf

@reference: Lorenz, Edward N. “Irregularity: A Fundamental Property of the 
            Atmosphere.” Tellus A 36A, no. 2 (1984): 98–110. 
            https://doi.org/10.1111/j.1600-0870.1984.tb00230.x.

@documentation:
    We make the distinction between two periods: the past and th future.
    The past is covered by one realization, the truth, and sampled by a limited
    number of observations
    The future is covered by many realizations, the ensemble, initialized
    from the present observation plus some noise.
    The dynamics obey the Lorenz 1984 system.
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
jM = 2

# Model parameters
a = 0.25
b = 4.0
F = 8.0
G = 1.0

# Time step for the numerical model expressed in days
dt = 1 / 6

# Run parameters
# When to start the truth, i.e. how many days before the initial time
day_past = -20

# Day of initialization
day_init  = 0.0 

# Day of verification
day_verif = 30#80

# Number of Monte-Carlo integrations
nMC     = 2
std_ini = 0.05 # standard deviation of Gaussian perturbation of initial state
std_obs =  0.1 # standard deviation of Gaussian perturbation for observations

# State label names
stl = ["$x_t$", "$y_t$", "$z_t$"]

# ==========================
# End of standard parameters
# ==========================
nday = day_verif - day_past
nt   = int(nday / dt )

# Creation of a time dimension
day = np.array([day_past + jt * dt for jt in range(nt + 1)])
# Time sequence
t   = np.arange(nt + 1)
t_init  = int((day_init  - day_past) / dt)
t_verif = int((day_verif - day_past) / dt)
# Index of truth, conveniently chosen as the first of the ensemble
jtruth = 0

# The system to be solved
#
# dX/dt = - Y^2 - Z^2 - aX + aF
# dY/dt =   XY  - bXZ - Y  + G
# dZ/dt =  bXY  + XZ  - Z

# Right-hand side of ODE as a function
# X is a column np array of size M
def f(t, M, X):
    # d is the day of year
    value = np.full(M, np.nan)
    
    value[0] = - X[1] ** 2 - X[2] ** 2 - a * X[0] + a *  F
    value[1] = X[0] * X[1] - b * X[0] * X[2] - X[1] + G * (1 - 2.0 * np.cos(2.0 * np.pi / (1.0 * 365.0) * t))
    value[2] = b * X[0] * X[1] + X[0] * X[2] - X[2]
    
    return value


# Initializing the state with "0"
# The first dimension is augmented by 1 to make space for the truth
X = np.full((1 + nMC, M, nt + 1), np.nan)
# Initializing the truth
X[jtruth, :, 0] = np.full(M, 0)

# Start integration
for jt in np.arange(1, nt + 1):
    # If we are before the initial time, let's just propagate the truth
    if day[jt] < day_init:
        X[jtruth, :, jt] = rk4vec(day[jt - 1], M, X[jtruth, :, jt - 1], dt, f)

       
    # Start Monte-Carlo integrations at initial time, sampled from truth
    elif day[jt] == day_init:
        
        # Update truth
        X[jtruth, :, jt] = rk4vec(day[jt - 1], M, X[jtruth, :, \
                 jt - 1], dt, f)
    
        # Generate initial condition from it
        Xini = X[jtruth, :, jt] 

        # Initialize all members                
        for jMC in np.arange(1, nMC + 1): # Ranging from 1 since 0 is the truth
            X[jMC, :, jt] = Xini + std_ini * np.random.randn(M)
            
    # If we are past the initial time, propagate from previous state
    else:
        for jMC in np.arange(nMC + 1):
            # Runge-Kutta 4th order method
            X[jMC, :, jt] = rk4vec(day[jt - 1], M, X[jMC, :, jt - 1], dt, f)


# Plot
# ----
fig, ax = plt.subplots(figsize = (8, 3), dpi = 300, \
                       constrained_layout = True)
plt.xlabel("Days")
plt.ylabel(stl[jM])
plt.xlim(day_past, 30)
plt.ylim(-3.0, 3.0)
plt.grid()

# Plot truth 
plt.plot(day, X[jtruth, jM, :], \
          lw = 3, \
          color = [0.5, 0.5, 0.5], label = "Solution de référence")
plt.legend()
plt.savefig("./fig0.png")

# Plot solution initialized from truth
plt.plot(day[t_init:], X[jtruth, jM, t_init:], \
          lw = 1, \
          color = [0.0, 1.0, 0.0], label = "Conditions initiales exactes")
plt.legend()
plt.savefig("./fig1.png")
# Plot solution initialized closed to truth
plt.plot(day[t_init:], X[1, jM, t_init:], \
          lw = 1, \
          color = [1.0, 0.0, 0.0], label = "Conditions initiales légèrement différentes")
plt.legend()
plt.savefig("./fig2.png")





















