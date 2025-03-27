# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import numpy as np
from utils import differencer
# Types and variables that will be plotted
types = {'Aerosol': ['PM25', 'AerMassNIT', 'AerMassNH4', 'AerMassPOA', 'AerMassBC', 'AerMassSO4'],
                  'O3' : ['SpeciesConc_O3'],
                  'NO2' : ['SpeciesConc_NO2']}
               
# --------------------------------------------------------------------------------------------------------------------
# Starting of the preprocessing, no need to modify anything after this

months = ['JAN', 'JUL']

# For loop to iterate through each type and through each variable withing each type
for type in types:
    for var in types[type]:
        fig, ax = plt.subplots(1, 2, figsize= [12, 4], subplot_kw={"projection": ccrs.EqualEarth(central_longitude=10)})
        plt.tight_layout()

        daSurf = [[],[]]
        for i, month in enumerate(months):

            daSurf[i] = differencer(type, month, var)

            if i==0:
                vmin = float(daSurf[i].values.min())
                vmax = float(daSurf[i].values.max())
            else:
                vmin = min(vmin, float(daSurf[i].values.min()))
                vmax = max(vmax, float(daSurf[i].values.max()))

        for i, month in enumerate(months):

            ax[i].add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
            ax[i].coastlines(resolution='50m', linewidth=0.5, color='black')

            # Plot the data for the current time step
            daSurf[i].plot(ax=ax[i], transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)

            # Add a title with the current time
            ax[i].set_title(f"{var}, monthly average {month}, level = {np.round(daSurf[i].lev.values, 3)}", fontsize=14)
        
       
        plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'monthly_average', 'pollutants_figures',f'{type}_{var}_timeaveraged.png'))
        
# Uncomment to display the plots on screen
# plt.show()