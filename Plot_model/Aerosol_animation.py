# sample_plot.py
#
# This script is an example showing how you can use [xarray] to access the
# data provided in the netCDF files and then use [matplotlib] and [cartopy]
# to plot a figure using the data
#
# For more information on these modules see:
#  https://docs.xarray.dev/en/stable/
#  https://matplotlib.org/index.html
#  https://scitools.org.uk/cartopy/docs/latest/
#
# — Flávio Quadros, 05/02/2020
#    - updated 07/03/2024

# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

# Open DataSet and print an overview of it
ds = xr.open_dataset(os.path.join('raw_data', 'model', 'Aerosol.JUL.OFF.nc4'))

# Select a DataArray
var = 'AerMassSO4'
da = ds[var]

# Set up the figure and axis
fig = plt.figure(figsize=[12, 6])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Function to update the plot for each time step
def update(frame):
    # Select the data for the current time step
    daSurf = da.isel(time=frame).isel(lev=0)
    
    # Clear the previous plot
    ax.clear()
    
    # Redraw the coastlines and borders for every frame
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
    ax.coastlines(resolution='50m', linewidth=0.5, color='white')
    
    # Plot the data for the current time step
    im = daSurf.plot(ax=ax, transform=ccrs.PlateCarree(), vmin=0, vmax=50, add_colorbar=False)
    
    # Add a title with the current time
    ax.set_title(f"{var} at time {da.time[frame].values}", fontsize=14)
    
    return [im]

# Create the animation
ani = FuncAnimation(fig, update, frames=len(da.time), interval=500, blit=False)

# Display the animation
plt.show()
