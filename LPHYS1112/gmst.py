# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:02:43 2019

@author: massonnetf
"""

import pandas as pd
import urllib
import matplotlib.pyplot as plt
import numpy as np
# Gaussian filter
def gsmooth(t_in, x_in, tau, t_out = None, method = None, check = False):
    """
    t_in      is a sorted N by 1 numpy array of independent coordinates
    x_in      is a N by 1 numpy array of data
    tau       is the decorrelation time scale for the smoothing window
              (in units of the time coordinate). When data is tau away from
              the center of the window, it is given a weight 1/e
    t_out  is a M by 1 numpy array of independent coordinates for smoothing. If
           no value is specified, a default number of 1000 points is taken
    method    is the method used to propagate data outside the window
    check     is a boolean to plot smoothing outside window
    
    returns
    
    gaussian smoothed version of x at t_out
    """
    
    # Convert to np array if necessary
    if type(t_in) is list:
        t_in = np.array(t_in)
    if type(x_in) is list:
        x_in = np.array(x_in)
    
    is_sorted = np.all(t_in[:-1] < t_in[1:])
    
    if not is_sorted:
        sys.exit("(gsmooth) t is not sorted")
        
    # Remove all NanS
    t_in = t_in[~np.isnan(x_in)]
    x_in = x_in[~np.isnan(x_in)]
   
    # 1. Mirror the data to sort the edges
    if method == "csym":  # Central symmetry

        x =  np.concatenate((np.array(2.0 * x_in[0] - x_in )[-1:0:-1], \
                             x_in                                    , \
                             np.array(2.0 * x_in[-1] - x_in)[-2::-1]     
                           ))
        t =  np.concatenate((np.array(2.0 * t_in[0] - t_in)[-1:0:-1],  \
                             t_in                                   ,
                             np.array(2.0 * t_in[-1] - t_in)[-2::-1]
                           ))
    elif method == None:
        pass
    else:
        sys.exit("(gsmooth) method unknown")
        
    if t_out == None:
        t_out = np.linspace(t[0], t[-1], 3* 1000) # 3 windows of 1000
    
    delta_t  = t_out[:, None] - t
    #weights  = np.exp(- np.abs(delta_t / tau)) # option exponential
    weights  = np.exp(- (delta_t / tau) ** 2)    # option gaussian
    weights /= np.sum(weights, axis = 1, keepdims=True)
    x_out    = np.dot(weights, x)

    start = np.where(t_out > t_in[0])[0][0]
    end   = np.where(t_out > t_in[-1])[0][0]

    if check:
        plt.figure("check")    
        plt.scatter(t_in,  x_in,  50, marker = "*", color = "black")
        plt.scatter(t, x, 10, marker = "*", color = "green")
        plt.plot(t_out, x_out, "red")
        plt.grid()
        plt.savefig("check.png", dpi = 300)
        plt.close("check")
    
    return t_out[start:end], x_out[start:end]

# Fetch the data
url = 'https://www.ncdc.noaa.gov/cag/global/time-series/globe/land_ocean/ytd/12/1880-2020.csv'  
filename = url.split("/")[-1]
urllib.request.urlretrieve(url, filename) 

# Read the data
df= pd.read_csv(filename, skiprows = 4, delimiter = ",")
years = df["Year"]
gmst  = df["Value"]

# Plot the data
fig1, ax1 = plt.subplots(1, 1, figsize = (5, 3), dpi = 300)
#plt.bar(years, gmst, width = 1)
for year, data in zip(years, gmst):
    color = plt.cm.inferno(int((data - min(gmst)) / (max(gmst) - min(gmst)) * 255.0))[:3]
    ax1.bar(year, data, width = 1, color = color )
    
ax1.grid()
ax1.plot((-1e9, 1e9), (0.0, 0.0), lw = 2, color = "k")
ax1.set_xlim(1870, 2022)
ax1.set_ylabel("$^\circ$ C")
ax1.set_title("Anomalies des températures moyennes\nglobales en surface (réf: 1901-2000)")
fig1.tight_layout()
    
fig1.savefig("./gmst1.png")


# Gaussian filter
years_out, gmst_smooth = gsmooth(np.array(years), np.array(gmst), tau = 2.0, method = "csym")
ax1.plot(years_out, gmst_smooth, linewidth = 3, color = [0.5, 0.5, 0.5])
fig1.savefig("./gmst2.png")

fig1.set_size_inches(4, 3)
ax1.set_xlim(1990, 2020)

fig1.tight_layout()
fig1.savefig("./gmst3.png")



# Sea ice 
url = "https://climexp.knmi.nl/data/iN_ice_extent.dat"
filename = url.split("/")[-1]
urllib.request.urlretrieve(url, filename) 
df= pd.read_csv(filename, skiprows = 7, delim_whitespace = True, header = None, engine = "python")
dfa = np.array(df.iloc[:, 1:])
dfa[dfa < 0.0] = np.nan
years2 = df[0][1:-1]
sie   = np.nanmean(dfa[1:-1, :], axis = 1)

fig2, ax2 = plt.subplots(1, 1, figsize = (5, 3))
for year, data in zip(years2, sie):
    color = plt.cm.RdYlBu(int((data - min(sie)) / (max(sie) - min(sie)) * 255.0))[:3]
    ax2.bar(year, data, width = 1, color = color )

ax2.set_title("Etendue annuelle moyenne de banquise arctique")
ax2.set_ylabel("million km$^2$")
ax2.set_ylim(9, 13)
ax2.grid()
fig2.tight_layout()
fig2.savefig("./sie1.png")



years_out, sie_smooth = gsmooth(np.array(years2), np.array(sie), tau = 2.0, method = "csym")
ax2.plot(years_out, sie_smooth, linewidth = 3, color = [0.5, 0.5, 0.5])
fig2.savefig("./sie2.png")
xlims = ax2.get_xlim()




fig1.savefig("./gmst4.png")


years = np.array(years)
gmst = np.array(gmst)
p = np.polyfit(np.arange(1992, 2001 + 1), gmst[(years >= 1992) * (years <= 2001)], 1)

ax1.set_xlim(xlims)
ax1.plot(np.arange(1992, 2001 + 1), np.arange(1992, 2001 + 1) * p[0] + p[1], "k--", lw = 3)
fig1.savefig("gmst5.png", dpi = 500)



p = np.polyfit(np.arange(2002, 2011 + 1), gmst[(years >= 2002) * (years <= 2011)], 1)
ax1.plot(np.arange(2002, 2011 + 1), np.arange(2002, 2011 + 1) * p[0] + p[1], "k--", lw = 3)
fig1.savefig("gmst6.png", dpi = 500)

