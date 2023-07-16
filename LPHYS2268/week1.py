from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
rundir = "/cofast/fmasson/PlaSim/plasim/CTRL/"

yearStart = 1
yearEnd   = 5223

years = np.arange(yearStart, yearEnd + 1)
nYears = len(years)

tasmean = np.full(nYears, np.nan)
for jYear, year in enumerate(years):
    print(year)
    f = Dataset(rundir + "/" + "CTRL." + str(year).zfill(5) + ".nc")

    tas = f.variables["tas"][:]
    lat = f.variables["lat"][:]
    lon = f.variables["lon"][:]

    dlam = 2 * np.pi / len(lon)
    dphi = 2 * np.pi / len(lat)

    nt = tas.shape[0]
    LON, LAT = np.meshgrid(lon, lat)

    tasmean[jYear] = np.mean([np.sum(np.squeeze(tas[jt, :, :]) * np.cos(LAT / 360.0 * 2 * np.pi)) / np.sum(np.cos(LAT / 360.0 * 2 * np.pi)) for jt in range(nt)])

    f.close()

fig, ax = plt.subplots(1, 1, figsize = (8, 4))
ax.plot(years, tasmean)
fig.tight_layout()
plt.savefig("./fig.png", dpi = 300)
