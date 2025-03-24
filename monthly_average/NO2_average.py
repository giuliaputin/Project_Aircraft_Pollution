# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import numpy as np

# Select variable you want to animate
# Choose: ['SpeciesConc_NO2']
var = 'SpeciesConc_NO2'

def differencer(type, month, var):
    dsoff = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." ,"raw_data", "model", f"{type}.{month}.OFF.nc4") )
    dson = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." ,"raw_data", "model", f"{type}.{month}.ON.nc4") )
    daoff = dsoff[var]
    daon = dson[var]

    daSurfoff = daoff.isel(lev=0).mean(dim = "time")
    daSurfon = daon.isel(lev=0).mean(dim = "time")
    daDiff = daSurfon - daSurfoff

    return daDiff
# --------------------------------------------------------------------------------------------------------------------
# Starting of the preprocessing, no need to modify anything after this

daSurf_JAN = differencer("NO2", "JAN", var)

daSurf_JUL = differencer("NO2", "JUL", var)

# Set up the figure and axis
fig, ax = plt.subplots(1, 2, figsize= [12, 6], subplot_kw={"projection": ccrs.EqualEarth(central_longitude=10)})
plt.tight_layout()

# ax[0] = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax[0].add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax[0].coastlines(resolution='50m', linewidth=0.5, color='black')

# Plot the data for the current time step
daSurf_JAN.plot(ax=ax[0], transform=ccrs.PlateCarree(), add_colorbar=False)

# Add a title with the current time
ax[0].set_title(f"{var}, monthly average JAN, level = {np.round(daSurf_JAN.lev.values, 3)}", fontsize=14)

# ax[1] = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax[1].add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax[1].coastlines(resolution='50m', linewidth=0.5, color='black')

# Plot the data for the current time step
daSurf_JUL.plot(ax=ax[1], transform=ccrs.PlateCarree(), add_colorbar=False)

# Add a title with the current time
ax[1].set_title(f"{var}, monthly average JUL, level = {np.round(daSurf_JUL.lev.values, 3)}", fontsize=14)

# Display the animation
plt.show()