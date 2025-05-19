# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import numpy as np

miscer = ["Met_LAI", "Met_PHIS"]
plotargs = [["Leaf Area Index", "BuGn", "leafareaindex"], ["Surface Geopotential Height", "terrain", "altitude"]]

# Which one would you like to plot: 0 Leaf Index, 1 Altitude
ploter = 1

# For loop to iterate through each type and through each variable withing each type
fig, ax = plt.subplots(1, 1, figsize= [9, 7], subplot_kw={"projection": ccrs.EqualEarth(central_longitude=10)})
plt.tight_layout()

ds = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." ,  "raw_data", "misc", "GEOSChem.StateMet.20190101_0000z.nc4") )
da = ds[miscer[ploter]]
damean = da.mean(dim = "time")
print(damean)

# daSurfoff = daoff.isel(lev=0).mean(dim = "time")
# daSurf[i] = daSurfoff


# vmin = float(da.values.min())
# vmax = float(da.values.max())

ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='black')

# Plot the data for the current time step
damean.plot(ax=ax,
            transform=ccrs.PlateCarree(),
            cmap=plotargs[ploter][1],
            cbar_kwargs={
                "label": plotargs[ploter][0],
                # "aspect": 40,
                "location": "bottom",
                "pad": 0.05,
                # "shrink": 0.4,
            })

# Add a title with the current time
# ax.set_title("Leaf Area Index", fontsize=14)


plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'backgroundgraphs', 'miscellaneous_figures', f'{plotargs[ploter][2]}.png'))

# Uncomment to display the plots on screen
# plt.show()