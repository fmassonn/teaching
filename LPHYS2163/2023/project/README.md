# LPHYS2163 Project #2
## How geostrophic is the real wind?

Here are the instructions for the second project of LPHYS2163.

### Objective 1. Producing the map of a variable

This first task is a technical one and ensures that you will have familiarized yourselves with basic tools. The preferred language for coding is python but alternatives are of course acccepted.

Your objective is to plot a map of the geopotential (in m$^2$/s$^2$) of the atmosphere at the 500 hPa level, on the 1st of December 2022 at 12:00 UTC. 

The data is already available in NetCDF format at `file:///Users/massonnetf/Dropbox/git/teaching/LPHYS2163/2023/project/data/sample_20221201_1200.nc`. This file is a "NetCDF" file, a common format used in geosciences. It features five variables:

- Temperature of the air
- Geopotential
- Zonal component of the wind
- Meridional component of the wind
- Relative vorticity

The variables are available at two pressure levels: 500 hPa (middle atmosphere) and 1000 hPa (surface)

You must first produce a map of the geopotential at the 500 hPa level like this one:


<p align="center">
<img src="./figs/z500_2022-12-01-1200.png">
</p>
