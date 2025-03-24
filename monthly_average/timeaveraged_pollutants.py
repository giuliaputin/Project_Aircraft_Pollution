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
types = {'Aerosol': ['PM25', 'AerMassNIT', 'AerMassNH4', 'AerMassPOA', 'AerMassBC'],
                  'O3' : ['SpeciesConc_O3'],
                  'NO2' : ['SpeciesConc_NO2']}


# It will now plot over: {'Aerosol': ['PM25', 'AerMassNIT', 'AerMassNH4', 'AerMassPOA', 'AerMassBC'],
#                  'O3' : ['SpeciesConc_O3'],
#                  'NO2' : ['SpeciesConc_O3']}                   

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

months = ['Jan', 'JUL']
# Set up the figure and axis


for type in types:
    for var in types[type]:
        fig, ax = plt.subplots(1, 2, figsize= [12, 4], subplot_kw={"projection": ccrs.EqualEarth(central_longitude=10)})
        plt.tight_layout()
    
        for i, month in enumerate(months):

            daSurf = differencer(type, month, var)
            ax[i].add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
            ax[i].coastlines(resolution='50m', linewidth=0.5, color='black')

            # Plot the data for the current time step
            daSurf.plot(ax=ax[i], transform=ccrs.PlateCarree(), add_colorbar=False)

            # Add a title with the current time
            ax[i].set_title(f"{var}, monthly average {month}, level = {np.round(daSurf.lev.values, 3)}", fontsize=14)
        
       


# Display the plots
plt.show()