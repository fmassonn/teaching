# to be executed on zenobe if the data is to be read
from netCDF4 import Dataset
import numpy as np
from seaice_commondiags import *
import matplotlib.pyplot as plt
import datetime
from datetime import timedelta
from matplotlib.dates import DateFormatter



root = "/storepelican/fmasson/APPLICATE/data/"

exp = "r"

loadRawData = False

nmembers = 2
members = np.arange(1, nmembers + 1)

myvars = ["siarean",]# "siareas", "siextentn", "siextents", "sivoln", "sivols"]
myvars = ["sivoln",]
nt = 731 # days
time = [datetime.date(2011, 1, 1) + timedelta(d) for d in range(nt)]

if loadRawData:
    
    # Read grid parameters
    f = Dataset("/storepelican/CLIMDATA/grid/mesh_mask_nemo.N3.6_ORCA1L75.nc", mode = "r")
    lat = f.variables["gphit"][0, :, :]
    lon = f.variables["glamt"][:][0, :, :]
    tmask= f.variables["tmaskutil"][:][0, :, :]
    e2t = f.variables["e2t"][0, :, :]
    e1t = f.variables["e1t"][0, :, :]
    f.close()


    siarean = list()
    siareas = list()
    siextentn = list()
    siextents = list()
    sivoln = list()
    sivols = list()
    
    for jmember, member in enumerate(members):
      file = root + exp + str(member).zfill(3) + "_20110101_20121231.nc"
      print(file)
      f = Dataset(file, mode = "r")
      sic = f.variables["siconc"][:] * 100.0 # in %
      siv = f.variables["sivolu"][:]
      f.close()
    
      # Compute statistics
      siarean.append(  compute_area(  sic, e1t * e2t ,                    mask = tmask * (lat > 0.0)))
      siareas.append(  compute_area(  sic, e1t * e2t ,                    mask = tmask * (lat < 0.0)))
      siextentn.append(compute_extent(sic, e1t * e2t , threshold = 15.0,  mask = tmask * (lat > 0.0)))
      siextents.append(compute_extent(sic, e1t * e2t , threshold = 15.0,  mask = tmask * (lat < 0.0)))
      sivoln.append(   compute_volume(siv, e1t * e2t ,                    mask = tmask * (lat > 0.0)))
      sivols.append(   compute_volume(siv, e1t * e2t ,                    mask = tmask * (lat < 0.0)))
    
    
    siarean = np.array(siarean).transpose()
    siareas = np.array(siareas).transpose()
    siextentn = np.array(siextentn).transpose()
    siextents = np.array(siextents).transpose()
    sivoln    = np.array(sivoln).transpose()
    sivols    = np.array(sivols).transpose()
    
    
    # Write to NetCDF
    f = Dataset("diags_" + exp  + "_20110101_20121231.nc", mode = "w", format = "NETCDF4_CLASSIC")
    
    time_dim = f.createDimension("time", None) # unlimited axis
    
    member_dim = f.createDimension("member", nmembers)
      
    out = f.createVariable("siarean",   np.float64,("time", "member",))
    out[:, :] = siarean
    out = f.createVariable("siareas",   np.float64,("time", "member",))
    out[:, :] = siareas
    out = f.createVariable("siextentn",   np.float64,("time", "member",))
    out[:, :] = siextentn
    out = f.createVariable("siextents",   np.float64,("time", "member",))
    out[:, :] = siextents
    out = f.createVariable("sivoln",   np.float64,("time", "member",))
    out[:, :] = sivoln
    out = f.createVariable("sivols",   np.float64,("time", "member",))
    out[:, :] = sivols
    
    f.close()
else:
    file = "./data/" + "diags_r_20110101_20121231.nc"
    f = Dataset(file, mode = "r")
    for var in myvars:
        vars()[var] = f.variables[var][:]
    f.close()


date_form = DateFormatter("%b %Y")


for var in myvars:
  series = vars()[var]
  fig, ax = plt.subplots(3, 1, figsize = (6, 9), dpi = 150)
  
  ax[0].plot(time, series, color = [0.5, 0.5, 0.5], lw = 0.1)
  ax[0].xaxis.set_major_formatter(date_form)
   
  ax[1].plot(time, np.std(series, axis = 1), color = "k")
  
  ax[2].plot(time, 1 / np.arange(1, nt + 1) * np.log(np.std(series, axis = 1) / np.std(series[0, :])), color = "k")
  for a in ax:
    for tick in a.get_xticklabels():
        tick.set_rotation(45)
  fig.tight_layout()
 

#)






