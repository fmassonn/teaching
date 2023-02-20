# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:02:43 2019

@author: massonnetf
"""

import pandas
import pandas as pd
import urllib
import matplotlib.pyplot as plt
import numpy as np

# Fetch the data
url = 'https://www.ncdc.noaa.gov/cag/global/time-series/globe/land_ocean/ytd/12/1880-2022.csv'  
filename = url.split("/")[-1]
urllib.request.urlretrieve(url, filename) 

# Read the data
df= pd.read_csv(filename, skiprows = 4, delimiter = ",")
years = df["Year"]
gmst  = df["Value"]

# Plot the data
fig1 = plt.figure("fig1", figsize = (5, 3), dpi = 300)
#plt.bar(years, gmst, width = 1)
for year, data in zip(years, gmst):
    color = plt.cm.inferno(int((data - min(gmst)) / (max(gmst) - min(gmst)) * 255.0))[:3]
    plt.bar(year, data, width = 1, color = color )
    
plt.grid()
plt.plot((-1e9, 1e9), (0.0, 0.0), lw = 2, color = "k")
plt.xlim(1970, 2024)
plt.ylabel("$^\circ$ C")
plt.title("Global mean surface temperature\nanomaly (ref: 1901-2000)")
plt.tight_layout()
    
plt.savefig("./gmst1.png")

# Linear fit, 1980-2010
istart = 90 # Year to start detrending
p = np.polyfit(years[istart:], gmst[istart:], 1)
fit = p[1] + p[0] * years[istart:]

plt.plot(years[istart:], fit, color = "cyan")
plt.savefig("./gmst2.png")

fig2 = plt.figure(figsize = (5, 3), dpi = 300)
ano = gmst[istart:] - fit

for j in np.arange(istart, len(years)):
    data = gmst[istart:] - fit
    color = plt.cm.RdBu_r(int((data[j] - min(data)) / (max(data) - min(data)) * 255.0))[:3]
    plt.bar(years[j], gmst[j] - fit[j], color = color)
plt.ylabel("°C")
plt.title("Anomaly of global mean surface temperature\n with respect to the linear trend line")
plt.grid()
plt.gca().set_axisbelow(True)
plt.tight_layout()
plt.savefig("./gmst3.png")
