# -*- coding: utf-8 -*-
"""
Created on Mon May 13 21:26:34 2019

@author: massonnetf
"""
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic
import numpy as np


f = Dataset("./sst.nc", mode = "r") 
lon = f.variables["longitude"][:]
lat = f.variables["latitude"][:]
sst = f.variables["sst"][:, :,:] - 273.15



map = Basemap(projection = "moll",  resolution = 'l', lon_0 = 0.0)
#map = Basemap(projection = "ortho", lat_0 = 0.0, lon_0 = 0.0)
# Map the lon, lat to the projection
lon, lat = np.meshgrid(lon, lat)
x, y = map(lon, lat)

clevs = np.arange(0.0, 30.0)
cmap = plt.cm.magma
mois = ["janvier", "février", "mars", "avril", "mai", "juin", "juillet", "août", "septembre", "octobre", "novembre", "décembre"]
for m in range(12):
    # Open figure
    plt.figure("fig", figsize = (6, 4))
    
    # Draw coastlines, country boundaries, fill continents, meridians & parallels
    map.drawcoastlines(linewidth = 0.25)
    map.fillcontinents(color = 'grey', lake_color = 'w')
    
    map.drawmeridians(np.arange(0, 360, 30))
    map.drawparallels(np.arange(-90, 90, 10))
    
    
    # Create a contourf object called "cs"
    cs = map.contourf(lon, lat, sst[m, :, :], clevs, cmap = cmap, \
                      latlon = True, extend = "both")
    cbar = map.colorbar(cs, location = 'bottom', pad = "5%")
    cbar.set_label("$^\circ$C")
    plt.title("Température moyenne de l'océan en\nsurface en " + mois[m] + " (1979-2018)")
    plt.tight_layout()
    #plt.savefig("sst1.png", dpi = 500)
    
    map.contour(lon, lat, sst[m, :, :], [26.5], colors = "k", linewidths = 2, latlon = True)
    plt.savefig("sst" + str(m + 1).zfill(2) + ".png", dpi = 500)
    plt.close("fig")
