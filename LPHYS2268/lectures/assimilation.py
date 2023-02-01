# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:59:35 2020

@author: massonnetf
"""
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
# Minimal example to illustrate multivariate data assimilation
# ------------------------------------------------------------

# Let u, v be the state vector, where u (v) is the zonal (meridional)
# component of the wind
#
# Let h(u, v) = sqrt(u^2+v^2) be the speed, that we can measure.
# Linearizing around a reference value u0, v0 gives
# H = 1 / h(u0, v0) [u0 v0]

# Reference state for linearization
u0 = 4
v0 = 4

H = 1 / np.sqrt(u0 ** 2 + v0 ** 2) * np.matrix([u0, v0])


# Observation of wind speed 
yo = 4.5
# Covariance matrix of wind speed
R  = (0.5) ** 2

# First guess state
xf = np.matrix([[3], [2]])
# Assumed variances
siguf = 0.7
sigvf = 0.3 
rhouv = 0.8
P = np.matrix([[siguf ** 2, rhouv * siguf * sigvf], 
             [rhouv * siguf * sigvf, sigvf ** 2]]) 


# Optimal gain
K = P * np.transpose(H) * (H * P * np.transpose(H) + R) ** (-1)


# Analysis
xa = xf + K * (yo - H * xf)
# Covariance
A = (np.matrix(np.eye(2)) - K * H) * P * (np.matrix(np.eye(2)) - 
     np.transpose(H) * np.transpose(K)) + K * R * np.transpose(K)


# Figures
# -------
fig, ax = plt.subplots(1, 2, figsize = (6, 3), dpi = 300)
# Dummy variable for plotting pdfs
xx = np.linspace(0.0, 10.0, 10000)

# Cosmetics
# Panel 1
ax[0].grid()
ax[0].set_title("Observation space")
ax[0].set_ylabel("PDF")
ax[0].set_xlim(0.0, 6.0)
ax[0].set_ylim(0.0, 1.5)
ax[0].set_xlabel("Wind speed [ms$^{-1}$]")
ax[0].set_ylabel("PDF [(ms$^{-1}$)$^{-1}$]")
# Panel 2
ax[1].grid()
ax[1].set_title("State (model) space")
ax[1].set_xlim(0.0, 6.0)
ax[1].set_ylim(0.0, 6.0)
ax[1].set_xlabel("Zonal wind u [ms$^{-1}$]")
ax[1].set_ylabel("Meridional wind v [ms$^{-1}$]")



# Panel 1 = observation space
# Plot obs and error PDF
ax[0].plot((yo, yo), (0, 1e9), "r--", label = "Observation + PDF")
pdf = 1 / np.sqrt(2 * np.pi * R) * np.exp(-0.5 * (xx - yo) ** 2 / R)
ax[0].plot(xx, pdf, "r")


ax[0].legend()
plt.tight_layout()
plt.savefig("./fig001.png")



# Panel 2 = state space


ax[1].scatter(xf[0, 0], xf[1, 0], 60, marker = "x", color = "grey",
  lw = 4, label = "First guess + PDFs")
# PDF for u
pdf = 1 / np.sqrt(2 * np.pi * siguf ** 2) * np.exp(-0.5 * (xx - xf[0, 0]) ** 2 / siguf ** 2)
ax[1].plot(xx, pdf, "grey")

# PDF for v
pdf = 1 / np.sqrt(2 * np.pi * sigvf ** 2) * np.exp(-0.5 * (xx - xf[1, 0]) ** 2 / sigvf ** 2)
ax[1].plot(pdf, xx, "grey")

ax[1].legend()
plt.savefig("./fig002.png")

# Plot small points for showing correlation
for i in range(1000):
    x = np.random.multivariate_normal(mean = np.array(xf.transpose())[0], cov = P)
    ax[1].scatter(x[0], x[1], 0.1, color = "lightgrey", zorder = 0)

ax[1].legend()
plt.savefig("./fig003.png")


# Project first guess to obs space
Hxf = (H * xf)[0, 0]
ax[0].plot((Hxf, Hxf), (0, 1e9), "k--", label = "First guess + PDF")
sigHef = np.sqrt(H * P * np.transpose(H))[0, 0]
pdf = 1 / np.sqrt(2 * np.pi * sigHef ** 2) * np.exp(-0.5 * (xx - Hxf) ** 2 / sigHef ** 2)
ax[0].plot(xx, pdf, "k")

ax[0].legend()
plt.savefig("./fig004.png")


# Compute analysis and its covariance matrix

# Plot analysis
# -------------
ax[1].scatter(xa[0, 0], xa[1, 0], 60, marker = "x", color = "green", lw = 4, label = "Analysis + PDFs")
pdf =  1 / np.sqrt(2 * np.pi * A[0, 0]) * np.exp(-0.5 * (xx - xa[0, 0]) ** 2 / A[0, 0])
ax[1].plot(xx, pdf, color = "green")

pdf =  1 / np.sqrt(2 * np.pi * A[1, 1]) * np.exp(-0.5 * (xx - xa[1, 0]) ** 2 / A[1, 1])
ax[1].plot(pdf, xx, color = "green")

# Plot small points for showing correlation
for i in range(1000):
    x = np.random.multivariate_normal(mean = np.array(xa.transpose())[0], cov = A)
    ax[1].scatter(x[0], x[1], 0.1, color = "lightgreen", zorder = 0, alpha = 0.5)

ax[1].legend()
plt.savefig("./fig005.png")

# In obs space
Hxa = (H * xa)[0, 0]
ax[0].plot((Hxa, Hxa), (0, 1e9), "g--", label = "Analysis + PDF")
sigHea = np.sqrt(H * A * np.transpose(H))[0, 0]
pdf =  1 / np.sqrt(2 * np.pi * sigHea ** 2) * np.exp(-0.5 * (xx - Hxa) ** 2 / sigHea ** 2)
ax[0].plot(xx, pdf, "g")


ax[0].legend()
plt.savefig("./fig006.png")