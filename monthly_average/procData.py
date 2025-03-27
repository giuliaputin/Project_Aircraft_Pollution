# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import numpy as np
from utils import *
from sklearn.preprocessing import StandardScaler
import time
start = time.time()
months = ['JAN', 'JUL']

types = {
    'Aerosol': {
        'PM25': ['nvPM'],
        'AerMassNIT': ['NO2'],
        'AerMassNH4': ['NO2'], # blimp BLOMP find something
        'AerMassPOA': ['nvPM', 'HC'],
        'AerMassBC': ['nvPM']
    },
    'O3': {
        'SpeciesConc_O3': ['NO2', 'HC', 'CO']
    },
    'NO2': {
        'SpeciesConc_NO2': ['NO2']
    }
}

emissions = xr.open_dataset(os.path.join(os.path.dirname(__file__), '..', 'raw_data', 'emissions', 'AvEmMasses.nc4'))
# emittants:
# NO2, HC, CO, nvPM
# Iterate through each pollutant
for type, vars in types.items():
    # Iterate through the emittants of each pollutant
    for var, emittants in vars.items():
        fig, ax = plt.subplots(1, 2, figsize= [12, 4], subplot_kw={"projection": ccrs.EqualEarth(central_longitude=10)})
        plt.tight_layout()
        
        vmin = float('inf')
        vmax = float('-inf')
        measures = []
        
        for month in months:
            pollutant = differencer(type, month, var)

            sum_emittants = adder(emissions, emittants).drop_sel(lat=68.5).drop_isel(lon=-1)
            # print(sum_emittants.sel(coords= sum_emittants.coords[sum_emittants.coords != (68.5, 48.12)]))
    
            print(pollutant)
            print(sum_emittants)
            sum_emittants_scaled = xr.DataArray(
                                    StandardScaler().fit_transform(sum_emittants.values),
                                    dims=sum_emittants.dims,
                                    coords=sum_emittants.coords,
                                    attrs=sum_emittants.attrs
                                )
            
            pollutant_scaled = xr.DataArray(
                StandardScaler().fit_transform(pollutant.values),
                dims=pollutant.dims,
                coords=pollutant.coords,
                attrs=pollutant.attrs
            )

            measure = pollutant / sum_emittants * 1e7
            print("\n\n MMMMM \n\n")
            print(measure)
            measures.append(measure)
            
            vmin = min(vmin, measure.values.min())
            vmax = max(vmax, measure.values.max())
            
        for i, month in enumerate(months):
            ax[i].add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
            ax[i].coastlines(resolution='50m', linewidth=0.5, color='black')

            # Plot the data for the current time step
            measures[i].plot(ax=ax[i], transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)

            # Add a title with the current time
            ax[i].set_title(f"{var}, monthly average {month}", fontsize=14)
            
end = time.time()
    
plt.show()



print(f"Process run in {end - start} s")