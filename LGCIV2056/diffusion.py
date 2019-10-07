# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 16:16:15 2019

@author: massonnetf
# Illustration of Fick's law with Brownian motion
"""

# Imports
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Cleaning up
plt.close("all")

# Simulation setup
N = 100000 # Number of particles
nt = 19 # Number of iterations
L = 3 # Half-length of the box
sig = 0.1 # standard deviation for Brownian motion

# Initialization: X and Y coordinates of points
X = np.full((N, nt), np.nan)
Y = np.full((N, nt), np.nan)

# Initial position
# Draw from a uniform distribution to constrain the particles to stay
# in the left part of the domain


x0  = - L * np.random.rand(N)
y0  =   L * (1 - 2 * np.random.rand(N))
        
X[:, 0] = x0
Y[:, 0] = y0
    
for jt in range(1, nt):
   print(jt)
   for jN in range(N):     
        # Iterate forward
        x_rand = sig * np.random.randn()
        y_rand = sig * np.random.randn()
        X[jN, jt] = X[jN, jt - 1] + x_rand if np.abs(X[jN, jt - 1] + \
              x_rand) < L else X[jN, jt - 1]
        Y[jN, jt] = Y[jN, jt - 1] + y_rand if np.abs(Y[jN, jt - 1] + \
              y_rand) < L else Y[jN, jt - 1]
    
   fig = plt.figure(dpi = 100)
   plt.plot((-L, L, L, -L, -L), (-L, -L, L, L, -L), "k-")
   plt.scatter(X[:, jt], Y[:, jt], 0.1)
   plt.title("t = " + str(jt).zfill(3))
   plt.xlim(-1.05 * L, 1.05 * L)
   plt.ylim(-1.05 * L, 1.05 * L)
   plt.savefig("./" + str(jt).zfill(4) + ".png")

   
   
#
## Draw from a normal distribution 
#
#X[:, 0] = sig * np.random.rand(N)
#Y[:, 0] = sig * np.random.rand(N)
#
#for jt in range(1, nt):
#    X[:, jt] = X[:, jt - 1] + sig * np.random.randn(N)
#    Y[:, jt] = Y[:, jt - 1] + sig * np.random.randn(N)
#    
#    fig = plt.figure(dpi = 100)
#
#    plt.scatter(X[:, jt], Y[:, jt], 0.01)
#    plt.title("t = " + str(jt).zfill(3))
#    plt.axis('equal')
#    plt.xlim(-4 ,4)
#    plt.ylim(-3, 3)
#
#    plt.savefig("./figs/diffusion_" + str(jt).zfill(4) + ".png")
#    plt.close()
# 
#    
