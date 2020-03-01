# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 15:38:31 2020

@author: massonnetf
"""

# Cox and Stephenson model
# explained in Hawkisn and sutton

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys

np.random.seed(0)
# Years
dt = 1
yearb, yeare = 1990, 2100
yinit = 2000
t = np.arange(yearb, yeare, dt) 

nt = len(t)


# Solving temperature increase for a simple model of the climate
# system with heat capacity.
#
# Model: C d(DT)/dt = DR [radiative imbalance] = DQ + lambda * DT
# with DQ assuming RCP-like forms, lambda the free parameter

# Climate system's heat capacity. From Boer, G. J., M. Stowasser, 
# and K. Hamilton. “Inferring Climate Sensitivity from Volcanic Events.” 
# Climate Dynamics 28, no. 5 (February 15, 2007): 481–502. 
# https://doi.org/10.1007/s00382-006-0193-x.

# Value from their Table 1. Average: 2.36e8 J / m2 / K
# =
# 7.48 W/m2/K year
C = 7.48

# Boltzmann constant, W/m2/K^4
sig = 5.67e-8

# Reference climate feedbak parameter
lambda_ref = -1.0

# Reference forcing (stabilization value)
dQ_ref = 4.5

# Percentiles (will be mirrored)
CIs = [75, 95, 97.5]
perc = [100 - c for c in CIs[-1::-1]] + CIs

# Open figure



# 
# ========================================
# ========================================
#

# General solution
# ----------------
dQ = dQ_ref  * (t - t[0] ) / (t[-1] - t[0])
fig0, ax0 = plt.subplots(1, 2, figsize = (8, 3), dpi = 300)

dT = 1 / C * np.exp(lambda_ref * t / C) * np.cumsum(dQ * np.exp(- lambda_ref /C * t) * dt)

ax0[0].plot(t, dQ, "green")
ax0[0].set_ylabel("Wm$^{-2}$")
ax0[1].plot(t, dT, "k")
ax0[0].grid()
ax0[0].set_title("Radiative forcing")
ax0[1].set_title("Temperature response\n($\lambda = $" 
                + str(lambda_ref) + " Wm$^{-2}K^{-1})$")
ax0[1].set_ylabel(r"$^\circ$ C")

ax0[1].grid()
plt.tight_layout()
plt.savefig("./series.png")

fig, ax = plt.subplots(3, 3, dpi = 300, figsize = (10, 7))
# Part 1: uncertainty due to model response
# -----------------------------------------

# Sigmoid-like shape to resemble RCP forcing
# Forcing
#dQ =  dQ_ref * 0.5 * (1 + np.tanh(0.01 * (t - (2300 + 1750) / 2))) # 0.5 * 12.5 * (1 + np.tanh(0.01 * (t - 2050)))
dQ = dQ_ref * (t - t[0] ) / (t[-1] - t[0])
# Sampling uncertainty in climate feedback parameter
nlambdas = 1000
lambdas = -1.0 + 0.2 * np.random.randn(nlambdas)

# Solution
dT = np.array([1 / C * np.exp(lam * t / C) * np.cumsum(dQ * np.exp(- lam /C * t) * dt) 
     for lam in lambdas])



dTperc = [np.percentile(dT, p, axis = 0) for p in perc]

# Lambda percentiles
lambdasperc = [np.percentile(lambdas, p) for p in perc]
ll = np.linspace(-2.0, 0.0, 1000)
kernel = stats.gaussian_kde(lambdas)
pdf = kernel(ll).T

ax[0, 0].plot(ll, pdf, color = "red")


for j, p in enumerate(perc[:int(len(perc)/2)]):
    color = [1.0,  0.7 ** (j + 1), 0.7 ** (j + 1)]

    ax[0, 1].fill_between(t, dTperc[j], dTperc[len(perc) - 1 -j], color = color, \
                     alpha = 0.7, lw = 0, zorder = 0 + j, label = str(int(100 - 2*p)) + " %" )

    ll = np.linspace(lambdasperc[j], lambdasperc[len(perc) - 1 -j])
    uu = kernel(ll).T
    ax[0, 0].fill_between(ll, np.zeros(len(ll)), uu, color = color, label = str(int(100 - 2*p)) + " %")


ax[0, 1].plot(t, np.mean(dT, axis = 0), "w-", lw = 2)
ax[0, 1].set_ylabel(r"$^\circ$ C")
ax[0, 1].set_title("Temperature response")
ax[0, 1].set_xlim(yinit - 10, yeare)
ax[0, 1].set_ylim(-0.5, 6.0)
ax[0, 1].legend()
ax[0, 1].grid()

ax[0, 0].set_title("Climate feedback parameter")
ax[0, 0].set_xlabel("$\lambda$ [Wm$^{-2}$K$^{-1}$]")
ax[0, 0].set_ylabel(r"Density [Km$^2$W$^{-1}$]")
ax[0, 0].grid()
ax[0, 0].legend()
plt.tight_layout()

# Fractional uncertainty
# Defined as 90 % CI / mean
FU1 =  100.0 * (np.percentile(dT, 95, axis = 0) - np.percentile(dT, 5, axis = 0) ) / np.mean(dT, axis = 0)
ax[0, 2].plot(t, FU1, color = "red", lw = 2)
ax[0, 2].set_ylabel("[%]")
ax[0, 2].grid()
ax[0, 2].set_title("Fractional uncertainty:\nmodel response")
ax[0, 2].set_xlim(yinit - 10, yeare)
ax[0, 2].set_ylim(0.0, 200.0)


# Part 2: forcing uncertainty
# ---------------------------
# linear + exponential uncertainty

nforcings = 1000
dQ_stoch = np.random.randn(nforcings)
dQs = [dQ_ref * (t - t[0] ) / (t[-1] - t[0]) +
        0.1 * dQS * (np.exp(3.0 * (t - yinit) / (t[-1] - yinit)) - 1) * (t > yinit)  for dQS in dQ_stoch]

# Solution
dT = np.array([1 / C * np.exp(lambda_ref * t / C) * np.cumsum(dQ * np.exp(- lambda_ref /C * t) * dt) 
     for dQ in dQs])

dTperc = [np.percentile(dT, p, axis = 0) for p in perc]

# Forcing percentiles
dQperc = [np.percentile(np.array(dQs), p, axis = 0) for p in perc]


for j, p in enumerate(perc[:int(len(perc)/2)]):
    color = [ 0.6 ** (j + 0.3), 0.9, 0.6 ** (j + 0.3)]

    ax[1, 1].fill_between(t, dTperc[j], dTperc[len(perc) - 1 -j], color = color, \
                     alpha = 0.6, lw = 0, zorder = 0 + j, label = str(int(100 - 2*p)) + " %" )

    ax[1, 0].fill_between(t, dQperc[j], dQperc[len(perc) - 1 - j], color = color, label = str(int(100 - 2*p)) + " %")
    


ax[1, 1].plot(t, np.mean(dT, axis = 0), "w-", lw = 2)
ax[1, 1].set_ylabel(r"$^\circ$ C")
ax[1, 1].set_title("Temperature response")
ax[1, 1].set_xlim(yinit - 10, yeare)
ax[1, 1].set_ylim(-0.5, 6.0)
ax[1, 1].legend()
ax[1, 1].grid()

ax[1, 0].set_title("Forcing")
ax[1, 0].grid()
ax[1, 0].legend()
plt.tight_layout()

# Fractional uncertainty
# Defined as 90 % CI / mean
FU2 = 100.0 * (np.percentile(dT, 95, axis = 0) - np.percentile(dT, 5, axis = 0) ) / np.mean(dT, axis = 0)
ax[1, 2].plot(t, FU2, color = "green", lw = 2)
ax[1, 2].set_ylabel("[%]")
ax[1, 2].grid()
ax[1, 2].set_title("Fractional uncertainty:\nscenario")
ax[1, 2].set_xlim(yinit - 10, yeare)
ax[1, 2].set_ylim(0.0, 200.0)


# Part 3: internal variability
# ----------------------------
nmembers = 1000
internalv = 0.2 * np.random.randn(nmembers, nt)

internalvperc = [np.percentile(np.array(internalv), p, axis = 0) for p in perc]

dQ =  dQ_ref * (t - t[0]) / (t[-1] - t[0])
dT = [1 / C * np.exp(lambda_ref * t / C) * np.cumsum(dQ * np.exp(- lambda_ref /C * t) * dt) + 
      internalv[jmembers, :] for jmembers in range(nmembers)]

dTperc = [np.percentile(dT, p, axis = 0) for p in perc]

ll = np.linspace(-1.0, 1.0, 1000)
kernel = stats.gaussian_kde(np.array(internalv).flatten())
pdf = kernel(ll).T
ax[2, 0].plot(ll, pdf, "blue")


for j, p in enumerate(perc[:int(len(perc)/2)]):
    color = [ 0.6 ** (j + 0.3), 0.6 ** (j + 0.3), 0.9]

    ax[2, 1].fill_between(t, dTperc[j], dTperc[len(perc) - 1 -j], color = color, \
                     alpha = 0.6, lw = 0, zorder = 0 + j, label = str(int(100 - 2*p)) + " %" )
    
    ll = np.linspace(np.mean(internalvperc[j]), np.mean(internalvperc[len(perc) - 1 -j]))
    uu = kernel(ll).T
    ax[2, 0].fill_between(ll, np.zeros(len(ll)), uu, color = color, 
      label = str(int(100 - 2*p)) + " %")

    
ax[2, 1].plot(t, np.mean(dT, axis = 0), "w-", lw = 2)
ax[2, 1].set_ylabel(r"$^\circ$ C")
ax[2, 1].set_title("Temperature response")
ax[2, 1].set_xlim(yinit - 10, yeare)
ax[2, 1].set_ylim(-0.5, 6.0)
ax[2, 1].legend()
ax[2, 1].grid()

ax[2, 0].set_xlim(-1.0, 1.0)
ax[2, 0].set_title("Internal variability")
ax[2, 0].set_xlabel(r"Temperature [" + "$^\circ$C]")
ax[2, 0].set_ylabel("Density [$^\circ$ C$^{-1}$]")
ax[2, 0].legend()

plt.tight_layout()


# Fractional uncertainty
# Defined as 90 % CI / mean
FU3 = 100.0 * (np.percentile(dT, 95, axis = 0) - np.percentile(dT, 5, axis = 0) ) / np.mean(dT, axis = 0)
ax[2, 2].plot(t, FU3, color = "blue", lw = 2)
ax[2, 2].set_ylabel("[%]")
ax[2, 2].grid()
ax[2, 2].set_title("Fractional uncertainty:\ninternal variability")
ax[2, 2].set_xlim(yinit - 10, yeare)
ax[2, 2].set_ylim(0.0, 200.0)


plt.savefig("./fig006.png")

plt.figure(figsize = (6, 3), dpi = 300)

plt.plot(t, FU1 ,"red", label = "model response")
plt.plot(t, FU2 ,"green", label = "scenario")
plt.plot(t, FU3 ,"blue", label = "internal variability")


plt.plot(t, FU1 + FU2 + FU3, "k", label = "total")
plt.xlim(yinit - 10, yeare)
plt.ylim(0.0, 250)
plt.grid()
plt.legend()
plt.savefig("./fig007.png")