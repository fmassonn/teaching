from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
rootdir = "/cofast/fmasson/PlaSim/plasim/PPP/"

nMemb = 5

yearStart = 3000
yearEnd   = 3009

years = np.arange(yearStart, yearEnd + 1)
nYears = len(years)
nDayPerYear = 360


dataPPP = np.full((nMemb, nDayPerYear * nYears), np.nan)


initialize=True

for jMemb in range(nMemb):
    print(jMemb)
    for jYear, year in enumerate(years):
        print(year)
        f = Dataset(rootdir + "PPP_" + str(jMemb + 1).zfill(3) + "/" + "PPP_" + str(jMemb + 1).zfill(3) + "." + str(year).zfill(5) + ".nc")


        tas = f.variables["tas"][:]

        lat = f.variables["lat"][:]
        lon = f.variables["lon"][:]

        LON, LAT = np.meshgrid(lon, lat)
        
        cosLat = np.cos(LAT / 360 * 2 * np.pi)

        tasMean = np.array([np.sum(tas[jt, :, :] * cosLat) / np.sum(cosLat) for jt in range(nDayPerYear) ])
        
        dataPPP[jMemb, jYear * nDayPerYear:(jYear + 1) * nDayPerYear] = tasMean

        f.close()

