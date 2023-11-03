import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime, timedelta

def plot_geopotential(time_index):
    # Open the NetCDF file
    with nc.Dataset('sample.nc') as ds:
        # Check if the provided time index is valid
        if time_index < 0 or time_index >= len(ds.dimensions['time']):
            raise ValueError("Invalid time index")
            
        # Find the index corresponding to 500 millibars
        level_index = np.where(ds['level'][:] == 500)[0]
        if len(level_index) == 0:
            raise ValueError("500 millibars level not found")
        level_index = level_index[0]
        
        # Read the geopotential data
        z = ds['z'][time_index, level_index, :, :]
        time_val = ds['time'][time_index]
        base_time = datetime.strptime(ds['time'].units, "hours since %Y-%m-%d %H:%M:%S.%f")
        time_stamp = base_time + timedelta(hours=int(time_val))
        print(time_stamp)
        # Apply scale factor and add offset
        #z = z * ds['z'].scale_factor + ds['z'].add_offset
        
        # Get latitude and longitude values
        latitudes = ds['latitude'][:]
        longitudes = ds['longitude'][:]
        
        # Create a meshgrid for latitudes and longitudes
        lon, lat = np.meshgrid(longitudes, latitudes)
        
        # Create the plot with Orthographic projection centered over Europe
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Orthographic(central_longitude=10, central_latitude=50)})
        
        # Set up map features
        ax.set_global()
        ax.coastlines(resolution='110m')
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.LAND, edgecolor='k')
 
        levels = np.arange(48000, 59000 + 100, 1000)       
        # Plot the geopotential data
        contour = ax.contourf(lon, lat, z, levels = levels, transform=ccrs.PlateCarree(), cmap='tab20b', extend='both')
        
        # Add color bar
        plt.colorbar(contour, label=f'Geopotential ({ds["z"].units})', orientation='horizontal', fraction=0.05, pad=0.15)
        
        # Add title
        plt.title(f'Geopotential at 500 millibars ('  + time_stamp.strftime("%d %b %Y  %H:%M") + ')')
	
        plt.tight_layout()        
        # Save figure
        yyyymmddhhmm = time_stamp.strftime("%Y-%m-%d-%H%M")
        plt.savefig('figs/z' + str(500) + "_" + yyyymmddhhmm + '.png')

        plt.close()

# Example usage
for i in np.arange(0, 24):
    plot_geopotential(i)  # Replace 0 with the desired time index

