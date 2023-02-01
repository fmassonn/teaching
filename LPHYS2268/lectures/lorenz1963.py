# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:10:05 2020

@author: massonnetf

@reference: Lorenz, Edward N. “Deterministic Nonperiodic Flow.” Journal of 
            the  tmospheric Sciences 20, no. 2 (March 1, 1963): 130–41. 
            https://doi.org/10.1175/1520-0469(1963)020<0130:DNF>2.0.CO;2.


"""

# Imports
# -------
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy import stats
from matplotlib.animation import PillowWriter
import matplotlib.animation as animation


plt.close("all")


# Numerical methods for solving the system 
# of ordinary differential equations (ODEs)
from numerical_methods import *

# Seed for reproducibility
np.random.seed(0)

# Model parameters

# State vector dimension
M = 3

# Parameter values
s = 10
b = 8 / 3
r = 28

# Forcing term a la Palmer
# Palmer, T. N. “Extended-Range Atmospheric Prediction and the Lorenz Model.”
# Bulletin of the American Meteorological Society 74, no. 1 (January 1, 1993): 
# 49–66. https://doi.org/10.1175/1520-0477(1993)074<0049:ERAPAT>2.0.CO;2.

forcing = 00.0

# Applies a projection to the coordinates to emulate a projection
# and thus a change in the shape of the attractor
proj = False

# Time stuff
# Time step for the numerical model expressed in days
dt = 1 / 200
# When to start the truth, i.e. how many days before the initial time
day_past = -20
# Day of initialization
day_init  = 0.0 

# Day of verification
day_verif = 0.5#1.5#0.5

# Number of Monte-Carlo integrations
nMC     = 500
std_ini = 0.35 # standard deviation of Gaussian perturbation of initial state

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

# Index of truth
jtruth = np.random.randint(nMC) 

force_Xini = True

# The system to be solved
#
# dX/dt = s (y -x)
# dY/dt = rx - y - xz
# dZ/dt =  xy - bz

# Right-hand side of ODE as a function
# X is a column np array of size M
def f(t, M, X):
        
    value = np.full(M, np.nan)
    
    value[0] = s * (X[1] - X[0])
    value[1] = r * X[0] - X[1] - X[0] * X[2] + forcing
    value[2] = X[0] * X[1] - b * X[2]
    
    return value


# The initial condition
X0 = np.full(M, 1) 

# The first dimension is augmented by 1 to make space for the truth
X = np.full((1 + nMC, M, nt + 1), np.nan)
# Initializing the truth
X[jtruth, :, 0] = np.full(M, np.random.randn())

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
        if force_Xini:
                Xini = np.array([ -7.2189591021895305 , -11.924935374229225,  
                                 14.585616894732182])
                Xini = np.array([-0.9131026,  0.7073493 , 19.05254639])
        else:
            Xini = X[jtruth, :, jt] 

        # Initialize all members                
        for jMC in np.arange(1, nMC + 1): # Ranging from 1 since 0 is the truth
            X[jMC, :, jt] = Xini + std_ini * np.random.randn(M)
            
    # If we are past the initial time, propagate from previous state
    else:
        for jMC in np.arange(nMC + 1):
            # Runge-Kutta 4th order method
            X[jMC, :, jt] = rk4vec(day[jt - 1], M, X[jMC, :, jt - 1], dt, f)


# Apply projection if needed
if proj:
    for jMC in range(nMC):
        for jt in range(nt + 1):
            #X[jMC, 0, jt] = 0.7 * X[jMC, 0, jt]
            theta = -10.0  / 360.0 * 2 * np.pi
            a, b = np.cos(theta) * X[jMC, 0, jt] - \
                            np.sin(theta) * X[jMC, 2, jt], \
                            np.sin(theta) * X[jMC, 0, jt] + \
                            np.cos(theta) * X[jMC, 2, jt] 
            X[jMC, 0, jt] = a
            X[jMC, 2, jt] = b
    colattr = [1.0, 0.5, 0.5]
else:
    colattr = [0.6, 0.6, 0.0]

fig, ax = plt.subplots(figsize = (6, 5), dpi = 100, \
                       constrained_layout = True)

ax = plt.axes(projection = '3d')
ax.view_init(30, -55)
ax.set_xlim(-20, 20)
ax.set_ylim(-30, 30)
ax.set_zlim(0, 50)

ax.plot3D(X[jtruth, 0, :], X[jtruth, 1, :], X[jtruth, 2, :],lw = 1, 
          alpha = 0.5, color = colattr)

plt.savefig("./fig1.png", transparent = False)

#[ax.scatter3D(X[jMC, 0, t_init:], X[jMC, 1, t_init:], 
#             X[jMC, 2, t_init:], "z", 0.5,
#             color = [0.0, 0.0, 0.0]  ) for jMC in range(1, nMC)]
[ax.plot3D(X[jMC, 0, t_init:], X[jMC, 1, t_init:], 
             X[jMC, 2, t_init:], lw = 0.5,
             color = [0.0, 0.0, 0.0], alpha = 0.05) for jMC in range(1, nMC)]

ax.scatter3D(X[jMC, 0, t_init], X[jMC, 1, t_init], X[jMC, 2, t_init], "z",
            150, color = [0.2, 0.2, 0.2], alpha = 0.8, lw = 5,  marker = "x")

plt.savefig("./fig2.png",transparent = False)




