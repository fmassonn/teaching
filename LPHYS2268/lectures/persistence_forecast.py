# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 15:38:50 2020

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

# Plots
# -----
fig, ax = plt.subplots(1, 1, figsize = (6, 3), dpi = 300)

# Cosmetics
ax.grid()
ax.set_axisbelow(True)
ax.set_xlim(1979.5, 2010.5)
ax.set_ylim(15.0, 20.0)
ax.set_ylabel("°C")
ax.set_title("Forecasts and verification of SST at Belgian coast")
ax.scatter(years, verif, marker = "s", color = [0.0, 0.0, 0.5], label = "Verification")
ax.legend()
plt.savefig("./fig000.png")

# Threshold to define event
thresh = 18.0
tt = ax.plot((-1e9, 1e9), (thresh, thresh), "k--", label = "Threshold")
ax.legend()
plt.savefig("./fig001.png")

aa = ax.fill_between((-1e9, 1e9, 1e9, -1e9), (thresh, thresh, 1e9, 1e9), 
                color = "green", alpha = 0.2, label = "Event", zorder = -1)
bb = ax.fill_between((-1e9, 1e9, 1e9, -1e9), (-1e9, -1e9, thresh, thresh), 
                color = "red", alpha = 0.2, label = "Non-event", zorder = -1)
ax.legend()
plt.savefig("./fig002.png")

ax.fill_between((-1e9, 1e9, 1e9, -1e9), (-1e9, -1e9, 1e9, 1e9), 
                color = "white", alpha = 1.0, zorder = 0)
aa.remove()
bb.remove()

aa = ax.scatter(years, forec, marker = "^", color = [0.8, 0.0, 0.0], label = "Forecast")
ax.set_xlim(2000.5, 2010.5)
ax.set_ylim(17, 20)
ax.legend()
plt.savefig("./fig003.png")


TT = np.linspace(10.0, 25.0, 10000)
plotleg = True
for year in years:
    pdf = normal(TT, mu = forec[year - years[0]], \
                          sigma = stdev[year - years[0]])
    if plotleg:
        leg = " Forecast PDF"
        plotleg = False
    else:
        leg = None
    mx = normal(forec[year - years[0]], mu = forec[year - years[0]], sigma = stdev[year - years[0]] )
   
    pdf1 = normal(TT[TT>forec[year - years[0]]], mu = forec[year - years[0]], \
                          sigma = stdev[year - years[0]])
    pdf2 = normal(TT[TT<forec[year - years[0]]], mu = forec[year - years[0]], \
                          sigma = stdev[year - years[0]])
    
    ax.fill_between(year + pdf1, np.full(len(pdf1), forec[year - years[0]]), 
                    TT[TT>forec[year - years[0]]],
                    color = [1.0, 0.0, 0.0], alpha = 0.2, lw = 0)
    ax.fill_between(year + pdf2, TT[TT<forec[year - years[0]]], 
                    np.full(len(pdf2), forec[year - years[0]]), 
                    color = [1.0, 0.0, 0.0], alpha = 0.2, lw = 0, label = leg)
aa.remove()
tt[0].remove()
ax.legend()
plt.savefig("./fig004.png")


# Re-add line
tt = ax.plot((-1e9, 1e9), (thresh, thresh), "k--", label = "Threshold")

ax.legend()
plt.savefig("./fig005.png")


# Scatter plot
fig, ax = plt.subplots(1, 1, figsize = (3, 3), dpi = 300)
ax.set_xlim(15.0, 20.0)
ax.set_ylim(15.0, 20.0)
ax.set_xlabel("Forecasted August SST [°C]")
ax.set_ylabel("Observed August SST [°C]")
ax.grid()
ax.plot((-1e9, 1e9), (-1e9, 1e9), "k--", label = "1:1 line")
ax.scatter(forec, verif, 20, marker = "s", color = [0.2, 0.2, 0.2], label = "Data")
correl = np.corrcoef(forec, verif)[0, 1]
MSE = np.mean((forec - verif) ** 2)
B = np.mean(forec - verif)



ax.legend()
plt.tight_layout()
plt.savefig("./fig100.png")


# Contingency table
# -----------------
n = len(years)
# YES - YES [Hits]
a = sum((forec >= thresh) * (verif >= thresh))
# YES - NO [False alarm]
b = sum((forec >= thresh) * (verif < thresh))
# NO - YES [Missed]
c = sum((forec <= thresh) * (verif >= thresh))
# NO - NO [Correct rejection]
d = sum((forec < thresh)  * (verif < thresh))

# Print performance indices
# -------------------------
PC = (a + d) / n
print("Proportion correct -->" + str(PC))

TS = a / (a + b + c)
print("Threat score -->" + str(TS))

OR = a * d / (b * c)
print("Odds ratio -->" + str(OR))

BIAS = (a + b) / (a + c)
print("Bias -->" + str(BIAS))


# Print it
fig, ax = plt.subplots(1, 1, figsize = (2, 2), dpi = 300)
ax.axis('off')
ax.set_xlim(-1, 2)
ax.set_ylim(-1, 2)
ax.text(0, 1.7, "YES",   ha = "center", va = "center")
ax.text(1, 1.7, "NO",    ha = "center", va = "center")
ax.text(0.5, 2.0, "Observed",   ha = "center", va = "center")
ax.text(-0.7, 1, "YES",  ha = "center", va = "center", rotation = 90)
ax.text(-0.7, 0, "NO",  ha = "center", va = "center", rotation = 90)
ax.text(-1, 0.5, "Forecast",   ha = "center", va = "center", rotation = 90)

ax.text(0, 1, str(a), ha = "center", va = "center")
ax.text(1, 1, str(b), ha = "center", va = "center")
ax.text(0, 0, str(c), ha = "center", va = "center")
ax.text(1, 0, str(d), ha = "center", va = "center")
[ax.plot((x, x), (-0.5, 1.5), "k-") for x in [-0.5, 0.5, 1.5]]
[ax.plot((-0.5, 1.5), (x, x), "k-") for x in [-0.5, 0.5, 1.5]]

plt.savefig("contingency.png")


# Probabilistic forecasts
# -----------------------

# If the forecast follows a normal centered on forec and with
# standard deviation stdev, then the probability that the variable exceeds
# some threshold is equal to the probability of (variable - forec) / stdev
# exceeds (thresh - forec) / stdev, the former beign N(0,1)
# This area is 1 - CDF 

prob = 1 - norm.cdf(((thresh - forec) / stdev))

fig, ax = plt.subplots(1, 1, figsize = (6, 3), dpi = 300)
# Cosmetics
ax.grid()
ax.set_axisbelow(True)
ax.set_xlim(1979.5, 2010.5)
ax.set_ylim(0.0, 100.0)
ax.set_ylabel("%")
ax.set_title("Probability of forecasted SST exceeding " + str(thresh) + " °C")
ax.bar(years, 100 * prob)
ax.set_xlim(2000.5, 2010.5)
ax.set_xticks(np.arange(2001, 2011))
plt.savefig("./fig006.png")

# Reliability diagram
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
for lb in np.arange(0, 1, step):
    ub = lb + step
    # Count how many observations of the event occurring occured
    # whenever it was forecast that the event would happen with
    # a probability between lb and lb + 0.1
    
    # How many forecasts fall in this class
    nf = sum((lb <= prob) * (prob < ub))
    
    rel = sum((lb <= prob) * (prob < ub) * (verif > thresh)) / nf

    ax.scatter(100 * (lb + ub) / 2, 100 * rel, 50, marker = "s", color = "k")
    
    if white:
        ax.fill_between((100 * lb, 100 * lb, 100 * ub), (0, 100, 100), 
                        color = "k", alpha = 0.2, lw = 0)
    white = not white
    
plt.savefig("./fig008.png")

# Brier score
verif_event = 1.0 * (verif > thresh)
BS = np.mean((prob - verif_event) ** 2)

p_clim = np.sum(verif_event) / len(years)
BS_ref = np.mean((p_clim - verif_event) ** 2)

def write_csv_SST():
    
    # ----------------
    # Data preparation
    # ----------------
    
    # Reading  raw data
    # -----------------
    # Full period (file)
    yearb, yeare = 1850, 2018
    file = "./data/cobe_sst_" + str(yearb) + "-" + str(yeare) + ".nc"
    # Extract grid-point SST data at the desired lon/lat coordinates
    f = Dataset(file, mode = "r")
    lat = f.variables["lat"][:]
    lon = f.variables["lon"][:]
    jy = np.where(lat == 51.5)[0][0]
    jx = np.where(lon == 3.5)[0][0]
    var = f.variables["sst"][:, jy, jx]
    f.close()
    
    
    # Creatoion of dataset
    
    # Time series of the verifying data
    # ---------------------------------
    
    # Years defining the period to train the forecast
    yearbf, yearef = 1900, 2015
    
    # Month predictand: August
    mon_tar  =  7 # Python convention
    
    # Time series
    verif    =  var[(yearbf - yearb) * 12 + mon_tar:(yearef - yearb) * 12 + 12:12]
    
    # Time series of the forecast data
    # Month predictor: June
    mon_pre = 5   # Python convention: 0 = January
    clim_len = 30 # Number of years for defining means 
    
    # Time axis
    years = np.arange(yearbf, yearef + 1)
    ny = len(years)
    
    # Forecast, consiting in a mean and a standard deviation for each  year
    # (Assuming gaussian)
    forec = [np.full(ny, np.nan), np.full(ny, np.nan)]
    
    for year in np.arange(yearbf, yearef + 1):
        print("")
        print("Forecasting year " + str(year))
        print("---------------------")
        
        # 1. Compute predictor climatology over past 30 years
        # ---------------------------------------------
        y2 = year - 1
        y1 = year - clim_len
        
        # First and last indices defining the period
        i1pre = (y1 - yearb) * 12 + mon_pre
        i2pre = (y2 - yearb) * 12 + mon_pre
        print("  Using climatology " + str(y1) + "-" + str(y2))
        
        pre_mean = np.mean(var[i1pre:i2pre + 1:12])
        print("  Predictor mean: " + str(pre_mean) + " °C")
        
        # 2. Compute this year's May anomaly
        # ----------------------------------
        i = (year - yearb) * 12 + mon_pre
        pre_ano = var[i] - pre_mean
        print("  Predictor anomaly for this year: " + str(pre_ano) + " °C")
        
        # 3. Compute target mean
        # ----------------------
        i1tar = (y1 - yearb) * 12 + mon_tar
        i2tar = (y2 - yearb) * 12 + mon_tar
        tar_mean = np.mean(var[i1tar:i2tar + 1:12])
        
        print("  Target mean: " + str(tar_mean) + " °C") 
        
        # 4. Compute uncertainty
        # ----------------------
        # Uncertainty comes from the fact that the means are subject to
        # sampling variability. If n years are available to estimate the
        # means, then the std on the mean is std / sqrt(n)
        # We add up variances, assuming May and september are not correlated
        
        std = np.sqrt(
              1 / len(var[i1pre:i2pre + 1:12]) * \
              np.var(var[i1pre:i2pre + 1:12]) + 
              1 / len(var[i1tar:i2tar + 1:12]) * \
              np.var(var[i1tar:i2tar + 1:12])
                     )
            
        forec[0][year - yearbf] = tar_mean + pre_ano
        forec[1][year - yearbf] = std
        
        print("  Prediction: " + str(forec[0][year - yearbf]) + " °C") 
        
    
    # Write to text file
    # ------------------
    with open('SST_NorthSea.csv', 'w', newline='\n') as f:
        writer = csv.writer(f)
        writer.writerow(["Year", "Observed " + months[mon_tar] + " SST", \
                         "Forecast " + months[mon_tar] + " SST (mean)", \
                         "Forecast " + months[mon_tar] + " SST (std)"])
        
        for year in np.arange(yearbf, yearef + 1):
            writer.writerow([year, \
                         "{:.3f}".format(verif[year - yearbf]),    \
                         "{:.3f}".format(forec[0][year - yearbf]), \
                         "{:.3f}".format(forec[1][year - yearbf])])
    
        f.seek(0, os.SEEK_END)
        f.seek(f.tell()-2, os.SEEK_SET)
        f.truncate()