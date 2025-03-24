# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import numpy as np

# From here you can set up the animation:
# Select either July (JUL) or January (JAN)
month = 'JUL'

# Select variable you want to animate
# Choose between: ['PM25', 'AerMassNIT', 'AerMassNH4', 'AerMassPOA', 'AerMassBC']
var = 'AerMassNIT'

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
# Open DataSet and print an overview of it

daSurf = differencer("Aerosol",month,var)

# Set up the figure and axis
fig = plt.figure(figsize=[12, 6])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='black')

# Plot the data for the current time step
daSurf.plot(ax=ax, transform=ccrs.PlateCarree(), add_colorbar=False)

# Add a title with the current time
ax.set_title(f"{var}, monthly average {month}, level = {np.round(daSurf.lev.values, 3)}", fontsize=14)

# Display the animation
plt.show()