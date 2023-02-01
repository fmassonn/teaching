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

Y = np.random.randn(M)
Z = 1.5 * Y + 1.0 * np.random.randn(M)



ax, fig = plt.subplots(1, 1, dpi = 300)
ax = plt.axes(projection = '3d')
ax.view_init(50, -45)

plt.tight_layout()
ax.plot3D(Y, Z, np.zeros(M), ".", markersize = 0.1, color = [0.3, 0.3, 0.3])
ax.set_zlim(0.0, 0.4)
ax.set_ylim(-5.0, 5)
ax.set_xlim(-5.0, 5)
ax.set_xlabel("$\mathbf{y}$")
ax.set_ylabel("$\mathbf{z}$")


# Fit density for Y
y_pdf = np.linspace(-8, 8, 1000)
kernelY = stats.gaussian_kde(Y)
pdf_Y = kernelY(y_pdf).T
ax.plot3D(y_pdf, - 4 * np.ones(len(y_pdf)), pdf_Y, color = [0.0, 0.0, 1.0], lw = 3,
          label = "$P(\mathbf{y})$")


# Fit density for Z
z_pdf = np.linspace(-8, 8, 1000)
kernelZ = stats.gaussian_kde(Z)
pdf_Z = kernelZ(z_pdf).T
ax.plot3D(4 * np.ones(len(z_pdf)), z_pdf, pdf_Z, color = [1, 0.5, 0.0], lw = 3,
          label = "$P(\mathbf{z})$")


# Conditional PDF of Y given Z in [0, 1]
ax.plot3D((-5.0, 5.0), (2.0, 2.0), (0.0, 0.0), color = [1, 0.5, 0.0], linestyle =  "--")
ax.plot3D((-5.0, 5.0), (3.0, 3.0), (0.0, 0.0), color = [1, 0.5, 0.0], linestyle =  "--")
YgZ = Y[(Z < 3.0) * (Z >= 2.0)]

# Fit density
kernelYgZ = stats.gaussian_kde(YgZ)
pdf_YgZ = kernel(y_pdf).T
ax.plot3D(y_pdf, - 4 * np.ones(len(y_pdf)), pdf_YgZ, 
          color = [0.0, 176/255, 240/255], lw = 3,
          label = "$P(\mathbf{y}|\mathbf{z})$")
ax.text3D(6.0, 2.0, 0.0, "$d\mathbf{z}$", color = [1, 0.5, 0.0])
x = [4, 4, 4, 4, 4]
y = [2, 3, 3, 2, 2]
z = [0,0, kernelZ(3.0).T, kernelZ(2.0).T, 0]
verts = [list(zip(x,y,z))]
#face = ax.add_collection3d(Poly3DCollection(verts, alpha = 0.5, facecolors = [1.0, 0.5, 0.0]))
pc = Poly3DCollection(verts, alpha = 0.5, facecolor = [1.0, 0.5, 0.0])
ax.add_collection3d(pc)

plt.legend()
plt.savefig("./fig003.png")


