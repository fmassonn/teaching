# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 09:47:13 2020

@author: massonnetf
"""
import numpy as np
import matplotlib.pyplot as plt

# Raw data
# X = km left as indicated by gauge (list of lists)
X = [[903, 883, 871, 864, 646, 593, 557, 447, 407, 333, 270, 201, 161,
99, 41], 
#[908, 872, 823, 577, 520, 470, 420, 233, 196, 129]
]
# Y = Counter
Y = [[63907, 63928, 63940, 63945, 64088, 64118, 64141, 64218, 64250,
64294, 64330, 64376, 64406, 64443, 64479],
#[63327, 63362, 63410, 63557, 63593, 63633, 63678, 63793, 63829, 63868]
]

# Convert into
# X: how much the counter said was left
# Y: how much was actually left, computed as the last mileage + 
#    what was said to be left (no better assumption) minus the current
x = np.array([item for sublist in X for item in sublist])

list_tmp = [[(Y[k][-1] + X[k][-1]) - Y[k][j] for j in range(len(Y[k]))] for k in range(len(Y))]
y =  np.array([item for sublist in list_tmp  for item in sublist])

fig, ax = plt.subplots(1, 1, figsize = (4, 4), dpi = 200)
ax.set_xlim(0, 1000)
ax.set_ylim(Y[0][0] - 100, Y[0][0] - 100 + 1000)
ax.set_xlabel("Predicted left mileage [km]")
ax.set_ylabel("Total mileage [km]")
ax.scatter(X, Y, 50, marker = "s", color = "k", label = "Data")
ax.grid()
ax.set_axisbelow(True)
ax.legend()
fig.tight_layout()
fig.savefig("scatter0.png")


fig, ax = plt.subplots(1, 1, figsize = (4, 4), dpi = 200)
ax.set_xlabel("Predicted left mileage [km]")
ax.set_ylabel("Actual left mileage [km]")
ax.grid()
ax.set_axisbelow(True)
ax.plot(xx, xx, "b--", label = "y=x")
ax.scatter(x, y, 50, marker = 's', color = "k", label = "Data")
ax.legend()
fig.tight_layout()
fig.savefig("scatter1.png")

# Fit line
xmean = np.mean(x)
ymean = np.mean(y)
xtild = x - np.mean(x)
ytild = y - np.mean(y)

a = np.sum(xtild * ytild) / np.sum(xtild ** 2)
b = ymean - a * xmean

xx = np.linspace(0, 1000)
ax.plot(xx, a * xx + b, "r-", label = "Regression")

# Plot theoretical line
# (y - y[0]) / (x - x[0]) = -1
yy = xx - x[0] + y[0]



plt.legend()
plt.savefig("scatter2.png")