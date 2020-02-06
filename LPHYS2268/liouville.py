# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 18:03:38 2020

@author: massonnetf
"""

# Imports
# -------
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy import stats
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Numerical methods for solving the system 
# of ordinary differential equations (ODEs)
from numerical_methods import *

np.random.seed(0)

# The assumed solution is X(t) = X(0) * exp(t/tau)
tau = 6.0


# X at time t
# Number of ensembles
M= 1000

fig, ax = plt.subplots(figsize = (6, 3), dpi = 300, \
                       constrained_layout = True)
ax = plt.axes(projection = '3d')
ax.view_init(40, 10)

tf = 5 # days
dt = 1 / 24 # time step
t = np.arange(0.0, tf, dt)
nt = len(t)
X0 = 2.0 + 0.5 * np.random.randn(M)

X = np.full((M, nt), np.nan)

# Plot one solution

ax.plot3D(t, X0[0] * np.exp(t / tau) , np.zeros(nt), color = [0.0, 176 / 255, 240 / 255],
              linewidth = 2, label = "One solution")
ax.set_xlabel("days")
ax.set_ylabel("x")
ax.set_zlim(0.0, 1.0)
ax.set_ylim(-0.0, 7)
plt.legend()
plt.tight_layout()
plt.savefig("./fig005a.png")

# Plot PDF of initial state
x_pdf = np.linspace(0, 7, 1000)
kernel0 = stats.gaussian_kde(X0)
pdf_X0 = kernel0(x_pdf).T
ax.plot3D(np.zeros(len(x_pdf)), x_pdf, pdf_X0, label = "$P(x_t, t)$")
plt.legend()
plt.savefig("./fig005b.png")


for jM in np.arange(M):
    if jM == 0:
        label = "All solutions compatible\nwith initial state"
    else:
        label = None
    X[jM, :] = X0[jM] * np.exp(t / tau)
    
    ax.plot3D(t, X[jM, :], np.zeros(nt), color = [0.0, 176 / 255, 240 / 255],
              linewidth = 0.1, label = label)
plt.legend()
ax.legend().get_lines()[-1].set_linewidth(1.0)
plt.tight_layout()
plt.savefig("./fig005c.png")

# PDF at later time
kerneldt = stats.gaussian_kde(X[:, -1])
pdf_Xdt = kerneldt(x_pdf).T
ax.plot3D(t[-1] * np.ones(len(x_pdf)), x_pdf, pdf_Xdt, label = "$P(x_{t+dt}, t+dt)$", color = [1.0, 0.5, 0.5])
plt.legend()
ax.legend().get_lines()[-2].set_linewidth(1.0)

plt.tight_layout()
plt.savefig("./fig005d.png")


# Plot conserved area under curve
y0, y1 = 2.0, 2.2
x0, x1 = 0.0, 0.0
z0, z1 = kernel0(y0).T, kernel0(y1).T

x = [x0, x0, x0, x0]
y = [y0, y1, y1, y0]
z = [0, 0,   z1, z0]

verts = [list(zip(x,y,z))]
pc = Poly3DCollection(verts, alpha = 0.5, facecolor = [0.0, 176 / 255, 240 / 255] )
ax.add_collection3d(pc)


# Plot end area under curve
ydt0, ydt1 = y0 * np.exp(t[-1]/tau), y1 * np.exp(t[-1]/tau) 
xdt0 = t[-1]
x = [xdt0, xdt0, xdt0, xdt0]
y = [ydt0, ydt1, ydt1, ydt0]
zdt0, zdt1 = kerneldt(ydt0).T, kerneldt(ydt1).T
z = [0, 0,   zdt1, zdt0]
verts = [list(zip(x,y,z))]
#face = ax.add_collection3d(Poly3DCollection(verts, alpha = 0.5, facecolors = [1.0, 0.5, 0.0]))
pc = Poly3DCollection(verts, alpha = 0.7, facecolor = [1.0, 0.5, 0.5] )
ax.add_collection3d(pc)

# Plot boundary trajectories
ax.plot3D(t, y0 * np.exp(t / tau), "w--", lw = 1)
ax.plot3D(t, y1 * np.exp(t / tau), "w--", lw = 1)

plt.savefig("./fig005e.png")


