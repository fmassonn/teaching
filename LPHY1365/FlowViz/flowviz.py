import numpy as np
import matplotlib.pyplot as plt

# This is flowviz, a python script that helps visualizing 2-D stationary (time-independent) flows


Lx = 2.0
Ly = 1.0

nx = 10 # Number of points in the x-direction
ny = 10 # Number of points in the x-direction

x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)

X, Y = np.meshgrid(x, y)

# --- User definitions ---

U = -3.0 * X

V = 2.0 * Y


# Plot quiver plot

fig, ax = plt.subplots(1, 1, dpi = 300)
ax.quiver(X, Y, U, V, color = "grey")
fig.savefig("./quiver.png")

# Seed particles
nParticles = 10

xParticles = np.full((nParticles, nTime), np.nan)
yParticles = np.full((nParticles, nTime), np.nan)

# Initialize particles




