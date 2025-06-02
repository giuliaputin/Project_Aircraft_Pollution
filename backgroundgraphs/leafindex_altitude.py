# Import modules
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import numpy as np

miscer = ["Met_LAI", "Met_PHIS"]
plotargs = [["Leaf Area Index", "BuGn", "leafareaindex"], ["Surface Geopotential Height [m]", "terrain", "altitude"]]

# Which one would you like to plot: 0 Leaf Index, 1 Altitude
ploter = 0

# For loop to iterate through each type and through each variable withing each type
fig, ax = plt.subplots(1, 1, figsize= [9, 7], subplot_kw={"projection": ccrs.EqualEarth(central_longitude=10)})
plt.tight_layout()

font = {'family' : 'DejaVu Sans',
        'size'   : 20}

matplotlib.rc('font', **font)

ds = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." ,  "raw_data", "misc", "GEOSChem.StateMet.20190101_0000z.nc4") )
da = ds[miscer[ploter]]
damean = da.mean(dim = "time")
print(damean)

# daSurfoff = daoff.isel(lev=0).mean(dim = "time")
# daSurf[i] = daSurfoff




ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='black')

# Plot the data for the current time step
# vmin = float(damean.values.min())
vmax = float(damean.values.max())

if ploter == 0:
    vmax = 7


damean.plot(ax=ax,
            transform=ccrs.PlateCarree(),
            cmap=plotargs[ploter][1],
            vmax = vmax,
            cbar_kwargs={
                "label": plotargs[ploter][0],
                # "aspect": 40,
                "location": "bottom",
                "pad": 0.01,
                # "shrink": 0.4,
            })

# Add a title with the current time
# ax.set_title("Leaf Area Index", fontsize=14)


plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'backgroundgraphs', 'miscellaneous_figures', f'{plotargs[ploter][2]}.png'), dpi=300)

# Uncomment to display the plots on screen
# plt.show()