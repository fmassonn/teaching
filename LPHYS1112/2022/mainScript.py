#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 17:35:56 2022

@author: massonnetf
"""
# Appearance
import matplotlib
matplotlib.rcParams['font.family'] = "Arial Narrow"

import matplotlib.ticker as mticker

# Libraries, imports, etc.
import numpy             as np
import matplotlib.pyplot as plt
import csv

from netCDF4 import Dataset
from scipy import interpolate 

# TO-DO
# Resolve issue units SSI


# Physical constants
h     = 6.62607015e-34 # [m 2 / kg  /   s]  Planck's constant
sigma = 5.670374419e-8 # [W   / m^2 / K^4]  Stefan-Boltzmann constant
c     = 299792458      # [m / s]            Speed of light
kB    = 1.380649e-23   # [J / K]            Boltzmann constant
dSunEarth = 1495000e6  # [m]                Distance Sun-Earth 
rEarth=6356e3          # [m]                Earth's radius
rSun  =696340e3        # [m]
b     =2.897771955e-3  # [m K]              Wien's displacement constant

# Functions
def PlanckLaw(lam, T):
    
    return (2 * h * c ** 2 / lam ** 5) * 1 / \
        (np.exp(h * c / (lam * kB * T)) - 1)

# Acronyms
# TOA: Top of Atmosphere
# OLR: Outgoing Longwave radiation

# Observed TOA OLR
OLR = 238.5 # [W/m2] in annual mean, global mean, 2000-2004
T_blackEarth = (OLR / sigma) ** (1/4)

# Read International Atmosphere Data
# Source: https://www.engineeringtoolbox.com/international-standard-atmosphere-d_985.html

fileCSV = open("./data/InternationalAtmosphere.csv")

csvReader = csv.reader(fileCSV, delimiter = ";")
[next(csvReader) for j in range(4)]
elevation = list()
vTemperature = list()
pressure = list()

for row in csvReader:
    #print(row)
    elevation.append(float(row[0]))
    vTemperature.append(float(row[1]))
    pressure.append(float(row[2]))

vTemperature = np.array(vTemperature)
elevation = np.array(elevation)
pressure = np.array(pressure)

fileCSV.close()




# Plot blackbody emission curve
# -----------------------------
Tsun  = 5800  
                           # Emission temperature
waveLength = np.linspace(0.1, 40, 10000) * 1e-6 # [m] wavelength as a variable 
                                            # (up to 6000 nm)

# Solid angle subtended by the Earth from the sun
# That's equal to the surface of a disc with Earth's diameter
# place on a sphere centered at the sun and with radius equal
# to Earth-Sun distance
angleSubtended = rEarth ** 2 / (dSunEarth ** 2) * 4 * np.pi

theoreticalSSI   = PlanckLaw(waveLength, Tsun) * 7.27e-5 / 1e9

jfig = 0
fig, ax = plt.subplots(1, 1, figsize = (4.5, 2.5))



# Plot observed solar irradiance (top of atmosphere)
fileIn = "./data/ssi_v02r01_monthly_s202101_e202112_c20220202.nc"
f = Dataset(fileIn, mode = "r")
actualSSI = np.mean(f.variables["SSI"][:], axis = 0)
waveLengthFile = f.variables["wavelength"][:]

a000 = ax.plot(waveLengthFile, actualSSI, label = "Observé, soleil")


ax.set_ylabel("Radiance spectrale\n[W m$^{-2}$ nm$^{-1}$]")
ax.set_xlabel("Longueur d'onde [nm]")
ax.set_xlim(0, 2000)
ax.legend()
ax.set_title("Spectre")
fig.tight_layout()

fig.savefig("./fig" + str(jfig).zfill(3) + ".png", dpi = 300) ; jfig += 1

a001 = ax.plot(waveLength * 1e9, theoreticalSSI, label = "Corps noir " + \
        str(Tsun) + " K", color = "orange")
ax.legend()
fig.savefig("./fig" + str(jfig).zfill(3) + ".png", dpi = 300) ; jfig += 1



# Wien's displacement law
# -----------------------
lambdaSun = b / Tsun

lambdaEarth = lambdaSun * Tsun / T_blackEarth 
theoreticalTSI = PlanckLaw(waveLength, T_blackEarth) * 7.27e-5 / 1e9
a002 = ax.plot(waveLength * 1e9, theoreticalTSI, label = "Corps noir " \
        + str(int(np.round(T_blackEarth))) + " K", color = "red")
ax.set_xlim(0, 18000)
ax.legend()
a000.pop(0).remove()
fig.savefig("./fig" + str(jfig).zfill(3) + ".png", dpi = 300) ; jfig += 1

# After normalisation
a001.pop(0).remove()
a002.pop(0).remove()

ax.plot(waveLength * 1e6, theoreticalSSI / np.max(theoreticalSSI), \
        color = "orange", label = "Corps noir " + \
        str(Tsun) + " K" )
ax.plot(waveLength * 1e6, theoreticalTSI / np.max(theoreticalTSI), \
        color = "red", label = "Corps noir " + \
        str(int(np.round(T_blackEarth))) + " K")
ax.set_xscale("log")
ax.set_xlabel("Longueur d'onde [$\mu$m]")
ax.set_ylabel("Radiance spectrale normalisée")
labels = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
ax.set_xticks(labels)
ax.set_xticklabels([str(l) for l in labels])
ax.set_ylim(0.0, 1.5)
ax.set_xlim(0.1, 40)
ax.legend(ncol = 2)
fig.savefig("./fig" + str(jfig).zfill(3) + ".png", dpi = 300) ; jfig += 1




# Read and plot absorptivity data
# -------------------------------



gases = ["H2O",      "CO2",     "CH4",     "O3",     "N2O"]
cols  = ["#00B5E6",  "#9F2D20", "#4A6D62", "#0035AD","#4E2029"]
fig, ax = plt.subplots(len(gases) + 1, 1, figsize = (5, 7))

tmpArray = np.full((len(gases), len(waveLength)), np.nan)

for jGas in range(len(gases)):
    wl = list()
    ab = list()
    fileCSV = open("./data/Absorptivity_" + gases[jGas] + ".csv")
    csvReader = csv.reader(fileCSV, delimiter = ";")
    next(csvReader); next(csvReader)
    for jrow, row in enumerate(csvReader):
        wl.append(float(row[0]) / 1e6) # Convert to m 
        ab.append(float(row[1]) * 100.0)  # Convert to %
        


    
    interp = interpolate.interp1d(np.array(wl), np.array(ab), kind = "linear", \
                                  bounds_error = False)
    interpAbsorpt = interp(waveLength)
    
    # Store that value for summming in the end
    tmpArray[jGas, :] = interpAbsorpt
    
    fileCSV.close()  


    #ax[jGas].fill_between(wl, ab, color = cols[jGas])
    ax[jGas].fill_between(waveLength * 1e6, interpAbsorpt, color = cols[jGas])
    #ax[jGas].plot(waveLength, interp(waveLength))

    
totalAbsorpt = np.sum(tmpArray, axis = 0)
ax[-1].fill_between(waveLength * 1e6, totalAbsorpt, color = "black")
  
for j, a in enumerate(ax[:-1]):
    a.set_xscale("log")
    a.set_xlim(0.1, 40)
    a.set_ylabel("%")
    #a.set_title(gases[j])
    a.text(0.5, 60, gases[j], color = cols[j], fontweight = "bold")
    a.set_ylim(0, 100)
    a.set_xticks(labels)
    a.set_xticklabels([str(l) for l in labels])

ax[-1].set_xscale("log")
ax[-1].set_xlim(0.1, 40)
ax[-1].set_ylabel("%")
ax[-1].text(0.5, 60, "Atmosphère", color = "black", fontweight = "bold")
ax[-1].set_ylim(0, 100)


fig.tight_layout()
fig.savefig("./fig" + str(jfig).zfill(3) + ".png", dpi = 300) ; jfig += 1


# Fig. with both spectra and absorptivities
fig, ax = plt.subplots(2, 1, figsize = (5, 3))
ax[0].plot(waveLength * 1e6, theoreticalSSI / np.max(theoreticalSSI), \
        color = "orange", label = "Corps noir " + \
        str(Tsun) + " K" )
ax[0].plot(waveLength * 1e6, theoreticalTSI / np.max(theoreticalTSI), \
        color = "red", label = "Corps noir " + \
        str(int(np.round(T_blackEarth))) + " K")
ax[0].set_xscale("log")

ax[0].set_ylabel("Radiance spectrale\nnormalisée")


ax[1].fill_between(waveLength * 1e6, totalAbsorpt, color = "black")
ax[1].set_xscale("log")
ax[1].set_xlim(0.1, 40)
ax[1].set_ylabel("Absorptivité [%]")
ax[1].text(0.5, 60, "Atmosphère", color = "black", fontweight = "bold")
ax[1].set_ylim(0, 100)
ax[1].set_xlabel("Longueur d'onde [$\mu$m]")
ax[1].set_xticks(labels)

fig.tight_layout()

fig.savefig("./fig" + str(jfig).zfill(3) + ".png", dpi = 300) ; jfig += 1


