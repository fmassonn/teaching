#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 14:24:52 2022

@author: massonnetf
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['font.family'] = "Arial Narrow"

p = np.linspace(0, 1, 1000)

pc = np.linspace(0, 1, 1000)

PC, P = np.meshgrid(p, pc)

def BS(X, Y):
    BS = Y + X ** 2 - 2 * X * Y 
    return BS

#BS = PC + P ** 2
#BS = PC + P


fig, _ = plt.subplots(2, 1, figsize = (4, 3), dpi = 300)

ax = plt.axes(projection = "3d")
ax.view_init(20, 250)

CS = ax.plot_surface(P, PC, BS(P, PC),  cmap = plt.cm.RdYlGn_r)

ax.set_xlabel("Forecast probability of event")
ax.set_ylabel("Baseline event frequency")
ax.set_zlabel("Brier Score")

ax.set_zlim(0.0, 1.5)
ax.clabel(CS, inline=1, fontsize=10)
cbar = fig.colorbar(CS)

fig.tight_layout(pad = 0.5)

xx = np.linspace(0, 1, 100)
ax.scatter3D(xx, xx, BS(xx, xx), s = 10, color = "magenta", marker = ".")

# baseline p = 0.30
ax.scatter3D(xx, np.ones(len(xx)), BS(xx, 0.30), color = "blue", marker = ".")



fig, ax= plt.subplots(1, 1, dpi = 300, figsize = (3, 2))
ax.plot(xx, xx ** 2, color = "red")
ax.plot(xx, (1 - xx) ** 2, color = "green")
ax.set_xlabel("Forecast probability")
ax.set_ylabel("Brier score")
fig.tight_layout()
fig.savefig("./figB.png")


fig, ax= plt.subplots(1, 1, dpi = 300, figsize = (3, 2))
ax.plot(xx, xx * (-np.log(xx)) - (1 - xx) * np.log(1-xx), color = "k")

ax.set_xlabel("Forecast probability")
ax.set_ylabel("Expected gain")
fig.tight_layout()
fig.savefig("./figLog.png")
