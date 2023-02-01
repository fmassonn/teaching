# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:41:00 2020

@author: massonnetf
"""
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

from ts_analyses import *



# Identification of a mode in SST


# Data
# (re-download  at ftp://ftp.cdc.noaa.gov/Datasets/COBE2/sst.mon.mean.nc)
file = "./data/cobe_sst_1850-2015.nc"

yearb = 1966
yeare = 2015
f = Dataset(file, mode = "r")
lat = f.variables["lat"][:]
lon = f.variables["lon"][:]
var = f.variables["sst"][(yearb-1850)*12:(yeare-1850)*12+12]
f.close()

lmask = 1 * var[0, :, :].mask

lon, lat = np.meshgrid(lon, lat)
nt, ny, nx = var.shape

lon0, lon1 = 120.0, 300.0
lat0, lat1 = -30.0, 30.0

fig = plt.figure(dpi = 300)
## Projection 
m = Basemap(projection='moll',lon_0= -180, lat_0 = 20, resolution='c')
#
#
x, y = m(lon, lat)
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))

levels = np.arange(0.0, 35.0)
c = m.contourf(x, y, var[-1, :, :],  latlon = False,
               cmap = plt.cm.inferno, levels = levels)
cb = m.colorbar(c, label='Sea surface temperature [$^\circ$ C]')
cb.set_ticks([np.arange(0.0, 40.0, 5)])
plt.tight_layout(pad = 50)
plt.savefig("./fig001.png")


# Draw box
m.plot((lon0, lon1), (lat0, lat0), \
       color = [0.0, 1.0, 1.0], linestyle = "--", lw = 3, latlon = True)
m.drawgreatcircle(lon1, lat0, lon1, lat1, color = [0.0, 1.0, 1.0], linestyle = "--", lw = 3)
m.plot((lon1, lon0), (lat1, lat1), \
       color = [0.0, 1.0, 1.0], linestyle = "--",  lw = 3, latlon = True)
m.drawgreatcircle(lon0, lat1, lon0, lat0, color = [0.0, 1.0, 1.0], linestyle = "--", lw = 3)
plt.savefig("./fig002.png")


# Decomposition
# -------------
# 1. Project physical data into matrix (the most annoying part)
# -------------------------------------------------------------
mask = 1 * (lon >= 120) * (lon <= 300.0) * (lat > -30.0) * \
        (lat <= 30.0) * (1 - lmask)

# State dimension
n_sta = np.sum(mask)
n_sam = nt 

# Create 2-D matrix
X_raw = np.empty((n_sta, n_sam))
X_raw[:] = np.nan


X_raw = var[:, mask == 1.0]


# 2. Pre-process the random vector, i.e. remove seasonality, trends, ...
# ----------------------------------------------------------------------
print("Pre-process")
X_raw_mean = np.mean(X_raw, axis = 0)
X = np.empty(X_raw.shape)

# Make JFM means and compute anomaly
 
Xjfm = np.array([np.mean(X_raw[jt:jt+3, :], axis = 0) for jt in np.arange(0, nt, 12)])
X = np.array([Xjfm[jt, :] - np.mean(Xjfm, axis = 0) for jt in range(Xjfm.shape[0])])

# 4. Do singular value decomposition
# ----------------------------------
print("SVD")
u, s, v = np.linalg.svd(X, full_matrices = 1, compute_uv = 1)
u = -u
v = -v
#  !!!! WATCH OUT !!!
#  this function returns u, s, v such that u * s * v = X (v not tranposed)

# 5. Analyze eigenvalues
# ----------------------
plt.figure()
plt.subplot(2, 1, 1)
plt.semilogy(range(len(s)), s ** 2, '.')
plt.subplot(2, 1, 2)
varian = s ** 2 / np.sum(s ** 2)
cumvar = np.cumsum(s ** 2) / np.sum(s ** 2)
plt.plot(range(len(s)), cumvar, '.')
plt.savefig("./eig.png")

# 6. Modes to retain
# ------------------
n_mod = np.max((0, np.sum((np.cumsum(s ** 2) / np.sum(s ** 2) < 0.9)))) # nb modes to retain))

# 7. Compute relevant time series ("factors")
# -------------------------------------------
#T = np.dot(u[:, :n_mod], np.diag(s[:n_mod]))
T = np.dot(u[:, :], np.diag(s[:]))
# 8. Map modes
# ------------
mode_out = np.empty((n_mod, ny , nx))
mode_out[:] = np.nan

for j_mod in np.arange(n_mod):
  out = np.empty((ny, nx))
  out[:] = np.nan
  out[mask == 1.0] = v[j_mod, :]
  mode_out[j_mod, :, :,] = out
  


for j_mod in np.arange(n_mod):
    fig = plt.figure(dpi = 300)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,120.,30.))
    m.drawmeridians(np.arange(0.,420.,60.))
    #levels = np.arange()
    c = m.contourf(x, y, mode_out[j_mod, :, :], latlon = False,
                   cmap = plt.cm.RdBu_r)
    plt.title("Mode " + str(j_mod + 1) + " (" + str(np.round(varian[j_mod]*100)) + " %)")
    cb = m.colorbar(c)
    plt.savefig("./fig" + str(j_mod).zfill(3) + "a.png")
    
    # Associated time series
    fig = plt.figure(dpi = 300, figsize = (5, 2))
    plt.plot(np.arange(yearb, yeare + 1 ), T[:, j_mod] / np.std(T[:, j_mod]), \
             color = [0.3, 0.3, 0.3], lw = 4)
    plt.grid()
    plt.tight_layout()
    plt.savefig("./fig" + str(j_mod).zfill(3) + "b.png")
    