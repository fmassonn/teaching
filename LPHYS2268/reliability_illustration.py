# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 09:02:24 2020

@author: massonnetf
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import scipy

np.random.seed(0)

# Illustration of reliability diagrams
# ------------------------------------

nt = 10000 # Number of re-forecasts

# Threshold. Above this value, we say that this event is realized
thresh =  1.5


# We take the obs as a realization from a PDF. This PDF represents
# all possible future states that could be compatible with the current
# state, given the level of precision we have in our models and obser-
# vations of initial state (and forcing). So this is really the
# "incompressible" part of the forecast.

# We take this PDF to vary from case to case, and be represented
# by a Gaussian distribution with varying mean and std

mean_obs = np.random.randn(nt)

# exponential to guarantee > 1
std_obs  = np.exp(0.5 * np.random.randn(nt))
    
# Draw obs from there
verif = mean_obs + std_obs * np.random.randn(nt)

# Binary variable to say if the event realized or not
verif_event = 1.0 * (verif > thresh)

# Then we have the forecasts. These are PDFs as well. 
# We have two types of forecast: those with a PDF that varies from one time 
# to the next and that resembles the observed; and climatological forecasts

# we allow the possibility to have a bias in the mean and in the std

# Climatological forecast or not?
forecast_type = "standard"

if forecast_type == "clim":
    mean_forecasts  = np.mean(verif) * np.ones(nt)
    std_forecasts   = np.std(verif)  * np.ones(nt)
elif forecast_type == "standard":
    bias    =   0.0   # Bias in the mean
    disp_fact = 2.0 # Dispersion factor: > 1.0 to inflate the spread of forecasts
                # (Make them overdispersive or underconfident)
                # 0 < disp_fact < 1.0 for underdispersive or overconfident
              
    mean_forecasts = (mean_obs + bias)
    std_forecasts = disp_fact * std_obs
elif forecast_type == "perfect deterministic":
    mean_forecasts = verif
    std_forecasts  = np.ones(nt) * 1e-3

# Event forecast probability
prob  = 1 - norm.cdf(((thresh - mean_forecasts) / std_forecasts))


# Compute reliability diagram
# ---------------------------
fig, ax = plt.subplots(1, 1, figsize = (4, 4), dpi = 300)
ax.grid()
ax.set_axisbelow(True)
ax.set_xlim(-5.0, 105.0)
ax.set_ylim(-5.0, 105.0)
ax.plot((0, 100), (0, 100), "k--", label = "1:1 line")
ax.set_xlabel("Forecast probability [%]")
ax.set_ylabel("Observed frequency [%]")
step = 0.1
white = True
step = 0.1

for lb in np.arange(0.0, 1.0, step):
    ub = lb + step
    # Count how many observations of the event occurring occured
    # whenever it was forecast that the event would happen with
    # a probability between lb and lb + 0.1
    
    # How many forecasts fall in this class
    nf = sum((lb <= prob) * (prob < ub))
    
    rel = sum((lb <= prob) * (prob < ub) * verif_event) / nf

    ax.scatter(100 * (lb + ub) / 2, 100 * rel, 50, marker = "s", color = "k")
    
    if white:
        ax.fill_between((100 * lb, 100 * lb, 100 * ub), (0, 100, 100), 
                        color = "k", alpha = 0.2, lw = 0)
    white = not white
    
    # Plot distribution of forecast probabilities
    dy = 8 * np.sign(rel - (lb + ub) / 2)
    #ax.text(100 * (lb + ub) / 2, 100 * rel + dy, str(int(100 * nf / nt)) + "%",
    #        ha = "center", va = "center")
    ax.bar(1 + 100 * (lb + ub) / 2 * 0.2, 100 * nf / nt * 0.2, bottom = 80,
           color = "k", width = 1.5)
    #ax.plot((0.0, 20, 20, 0.0, 0.0), (80, 80, 100, 100, 80), "k-", lw = 1)

plt.tight_layout()
plt.savefig("./reliability.png")

# Plot fraction of observed frequencies following each class of event
fig, ax = plt.subplots(1, 1, figsize = (4, 4), dpi = 100)
ax.grid()
ax.set_axisbelow(True)
ax.set_xlim(-5.0, 105.0)

ax.hist(100 * prob[verif_event == 1], bins = 100 * np.arange(0.0, 1.0 +
        step, step), 
        density = True, alpha = 0.5, label = "When the event happened")
ax.hist(100 * prob[verif_event == 0], bins = 100 * np.arange(0.0, 1.0 + 
        step, step), 
        density = True, alpha = 0.5, label = "When the event did not happen")
ax.set_title("Distribution of forecasted probabilities")
ax.set_xlabel("Forecast probability [%]")
ax.set_ylabel("PDF")
ax.legend()


# Plot one as example
# -------------------
xx = np.linspace(-6.0, 6.0, 10000)

fig, axs = plt.subplots(1, 3, figsize = (12, 4  ), dpi = 300)
for j in range(3):
    ax = axs[j]
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_ylabel("PDF")
    ax.set_ylim(0.0, 1.0)
    ax.set_xlim(-2.0, 4.0)
    ax.set_xlabel("Physical variable")
    ax.set_title("Forecast #" + str(j + 1) + " (out of " + str(nt) + ")" + 
                 "\nBias = " + str(bias) + ", dispersion factor = " + 
                 str(disp_fact))
    
    # Plot OBS PDF
    pdfobs = 1.0 / np.sqrt(2 * np.pi * std_obs[j] ** 2) * \
           np.exp(- (xx - mean_obs[j]) ** 2 / (2 * std_obs[j] ** 2))
           
    ax.plot(xx, pdfobs, 
                    color = [0.0, 0.8, 0.8], lw = 4, ls = "--",
                    label = "PDF generating the observation")
   
    
    # Plots verifying obs
    ax.plot((verif[j], verif[j]), (-1e9, 1e9), color = [0.0, 0.8, 0.8],
            lw = 3, 
            label = "Verifying observation", ls = ":")
    
    if j == 0:
        ax.legend()
    plt.tight_layout()
    plt.savefig("./fig000_" + str(j).zfill(3) + ".png")
    

    
    # Plot forecast PDF
    pdf = 1.0 / np.sqrt(2 * np.pi * std_forecasts[j] ** 2) * \
           np.exp(- (xx - mean_forecasts[j]) ** 2 / 
                  (2 * std_forecasts[j] ** 2))
    ax.plot(xx, pdf, color = [0, 0, 0], label = "Forecast PDF", lw = 7, 
            zorder = 1)

    if j == 0:
        ax.legend()
    plt.savefig("./fig001_" + str(j).zfill(3) + ".png")
    
    ax.plot((thresh, thresh), (-1e9, 1e9), color = [0.0, 0.0, 0.0], 
            ls = "-",
            label = "Threshold defining the event", lw = 3)
    ax.plot((0, 0), (-1e9, 1e9), color = "k", lw = 1)
    
    if j == 0:
        ax.legend()
    plt.savefig("./fig002_" + str(j).zfill(3) + ".png")
    
    
    # Color to the right of the threshold
    ax.fill_between(xx[xx >  thresh], pdf[xx > thresh], 
                    color = [0.0, 0.5, 0.0], alpha = 0.7, 
                    label = "Forecasted event probability")
    
    # Color to the left of the threshold
    ax.fill_between(xx[xx <=  thresh], pdf[xx <= thresh], 
                    color = [0.5, 0.0, 0.0], alpha = 0.6, 
                    label = "Forecasted non-event probability")
    
    ax.text(2.9, 0.2, str(int(100.0 * prob[j])) + \
            "%", color = [0.0, 0.5, 0.0], fontsize = 16, fontweight = "bold")

    if j == 0:
        ax.legend()
    plt.savefig("./fig003_" + str(j).zfill(3) + ".png")
    
    # Plot background as green or red to mean it happened or not
    if verif_event[j] == 1:
        color = "green"
    else:
        color = "red"
        
    ax.fill_between((-10, 10, 10, -10, -10), (0, 0, 1, 1, 0), color = color,
                    alpha = 0.1, zorder = -1000)
    
    if j == 0:
        ax.legend()
    plt.savefig("./fig004_" + str(j).zfill(3) + ".png")
    
    # Plot PDF from which obs is drawn
    if j == 0:
        ax.legend()

plt.savefig("./fig.png")