# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 22:46:39 2020

@author: massonnetf
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import calendar
import csv
import os

from scipy.stats import norm

# Illustration of four types of forecasts
def normal(xx, mu = 0, sigma = 1):
    return 1.0 / np.sqrt(2 * np.pi * sigma ** 2) * np.exp(- ((xx - mu)  \
                / (np.sqrt(2) * sigma)) ** 2)

# Month labels
months = [calendar.month_name[i + 1] for i in range(12)]

# Sample data set for testing various concepts in forecast verification
# and calibration

# F. Massonnet, March 2020.

# Read input file
# ---------------
years  = list()
verif = list()
forec = list()
stdev = list()

file = "./data/SST_NorthSea.csv"
with open(file, 'r') as f:
    reader = csv.reader(f, delimiter = ",")
    # Skip first row with labels
    next(reader)
    for row in reader:
        years.append(int(  row[0]))
        verif.append(float(row[1]))
        forec.append(float(row[2]))
        stdev.append(float(row[3]))
        
years = np.array(years)                   
verif = np.array(verif)
forec = np.array(forec)
stdev = np.array(stdev)


# 1. Biased low forecast
forec1 = forec - 2.0
forec2 = 0.2 * (forec - np.mean(forec)) + np.mean(forec)
forec3 = forec - np.polyval(np.polyfit(years, forec, 1), years) + np.mean(forec)
# Plots
# -----
fig, ax = plt.subplots(1, 1, figsize = (6, 3), dpi = 300)

# Cosmetics
ax.grid()
ax.set_axisbelow(True)
#ax.set_xlim(1989.5, 2010.5)
ax.set_ylim(10.0, 20.0)
ax.set_ylabel("°C")
ax.set_title("Forecasts and verification of SST at Belgian coast")
ax.scatter(years, verif, marker = "s", color = [0.0, 0.0, 0.5], 
           label = "Verification")
ax.scatter(years, forec1, marker = "^", color = [0.8, 0.0, 0.0],
           label = "Forecast")
ax.legend()
plt.savefig("./fig000.png")

fig, ax = plt.subplots(1, 1, figsize = (4, 4), dpi = 300)
ax.scatter(forec1, verif, marker = ".", color = "k")
ax.plot((0.0, 100), (0.0, 100), "b--", label = "y=x")
ax.grid()
ax.set_xlim(13.0, 20.0)
ax.set_ylim(ax.get_xlim())
ax.set_ylabel("Verification SST [°C]")
ax.set_xlabel("Forecast SST [°C]")
ax.set_title("Bias in the mean")
ax.legend()
plt.tight_layout()
plt.savefig("./fig001.png")

fig, ax = plt.subplots(1, 1, figsize = (4, 4), dpi = 300)
ax.scatter(forec1, verif, marker = ".", color = [0.9, 0.9, 0.9], 
           label = "Raw forecasts")
ax.scatter(forec1 + np.mean(verif - forec1), verif, marker = ".", 
           color = "orange", label = "Mean bias-corrected forecasts")
ax.plot((0.0, 100), (0.0, 100), "b--", label = "y=x")
ax.grid()
ax.set_xlim(13.0, 20.0)
ax.set_ylim(ax.get_xlim())
ax.set_ylabel("Verification SST [°C]")
ax.set_xlabel("Forecast SST [°C]")
ax.set_title("Bias in the mean removed")
ax.legend()
plt.tight_layout()
plt.savefig("./fig002.png")

# -------------


# Plots
# -----
fig, ax = plt.subplots(1, 1, figsize = (6, 3), dpi = 300)

# Cosmetics
ax.grid()
ax.set_axisbelow(True)
#ax.set_xlim(1989.5, 2010.5)
#ax.set_ylim(14.0, 20.0)
ax.set_ylabel("°C")
ax.set_title("Forecasts and verification of SST at Belgian coast")
ax.scatter(years, verif, marker = "s", color = [0.0, 0.0, 0.5], 
           label = "Verification")
ax.scatter(years, forec2, marker = "^", color = [0.8, 0.0, 0.0],
           label = "Forecast")
ax.legend()
plt.savefig("./fig003.png")

fig, ax = plt.subplots(1, 1, figsize = (4, 4), dpi = 300)
ax.scatter(forec2, verif, marker = ".", color = "k")
ax.plot((0.0, 100), (0.0, 100), "b--", label = "y=x")
ax.grid()
ax.set_xlim(13.0, 20.0)
ax.set_ylim(ax.get_xlim())
ax.set_ylabel("Verification SST [°C]")
ax.set_xlabel("Forecast SST [°C]")
ax.set_title("Bias in the variability")
ax.legend()
plt.tight_layout()
plt.savefig("./fig004.png")

fig, ax = plt.subplots(1, 1, figsize = (4, 4), dpi = 300)
ax.scatter(forec2, verif, marker = ".", color = [0.9, 0.9, 0.9], 
           label = "Raw forecasts")
ax.scatter((forec2 - np.mean(forec2)) / np.std(forec2) * np.std(verif) + 
           np.mean(forec2), verif, marker = ".", 
           color = "orange", label = "Mean bias-corrected forecasts")
ax.plot((0.0, 100), (0.0, 100), "b--", label = "y=x")
ax.grid()
ax.set_xlim(13.0, 20.0)
ax.set_ylim(ax.get_xlim())
ax.set_ylabel("Verification SST [°C]")
ax.set_xlabel("Forecast SST [°C]")
ax.set_title("Bias in the variability removed")
ax.legend()
plt.tight_layout()
plt.savefig("./fig005.png")



# Plots
# -----
fig, ax = plt.subplots(1, 1, figsize = (6, 3), dpi = 300)

# Cosmetics
ax.grid()
ax.set_axisbelow(True)
#ax.set_xlim(1989.5, 2010.5)
#ax.set_ylim(14.0, 20.0)
ax.set_ylabel("°C")
ax.set_title("Forecasts and verification of SST at Belgian coast")
ax.scatter(years, verif, marker = "s", color = [0.0, 0.0, 0.5], 
           label = "Verification", alpha = 0.5)
ax.scatter(years, forec3, marker = "^", color = [0.8, 0.0, 0.0],
           label = "Forecast", alpha = 0.5)
ax.plot(years, np.polyval(np.polyfit(years, verif, 1), years), color = [0.0, 0.0, 0.5],
        label = "Trend", lw = 2)
ax.plot(years, np.polyval(np.polyfit(years, forec3, 1), years), color = [0.8, 0.0, 0.0],
        label = "Trend", lw = 2)
ax.legend(ncol = 2)
plt.savefig("./fig006.png")

fig, ax = plt.subplots(1, 1, figsize = (4, 4), dpi = 300)
ax.scatter(forec3, verif, marker = ".", color = "k")
ax.plot((0.0, 100), (0.0, 100), "b--", label = "y=x")
ax.grid()
ax.set_xlim(13.0, 20.0)
ax.set_ylim(ax.get_xlim())
ax.set_ylabel("Verification SST [°C]")
ax.set_xlabel("Forecast SST [°C]")
ax.set_title("Bias in the trend")
ax.legend()
plt.tight_layout()
plt.savefig("./fig007.png")

#fig, ax = plt.subplots(1, 1, figsize = (4, 4), dpi = 300)
#ax.scatter(forec2, verif, marker = ".", color = [0.9, 0.9, 0.9], 
#           label = "Raw forecasts")
#ax.scatter((forec2 - np.mean(forec2)) / np.std(forec2) * np.std(verif) + 
#           np.mean(forec2), verif, marker = ".", 
#           color = "orange", label = "Mean bias-corrected forecasts")
#ax.plot((0.0, 100), (0.0, 100), "b--", label = "y=x")
#ax.grid()
#ax.set_xlim(13.0, 20.0)
#ax.set_ylim(ax.get_xlim())
#ax.set_ylabel("Verification SST [°C]")
#ax.set_xlabel("Forecast SST [°C]")
#ax.set_title("Bias in the variability removed")
#ax.legend()
#plt.tight_layout()
#plt.savefig("./fig005.png")