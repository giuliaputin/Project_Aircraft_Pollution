import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import numpy as np
from utils import *
from sklearn.preprocessing import StandardScaler
import time
from cartopy.io.shapereader import natural_earth
import shapely.geometry as sgeom
import cartopy.io.shapereader as shpreader
from shapely.prepared import prep

start = time.time()
months = ['JAN', 'JUL']

types = {
    'Aerosol': {
        'PM25': ['nvPM'],
        'AerMassNIT': ['NO2'],
        'AerMassNH4': ['NO2'],
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

for type_, vars in types.items():
    for var, emittants in vars.items():
        fig, ax = plt.subplots(1, 2, figsize=[12, 3], subplot_kw={"projection": ccrs.PlateCarree()})
        plt.tight_layout()
        pollutants_lst = []
        emittants_lst = []
        
        for month in months:
            pollutant = differencer(type_, month, var)
            sum_emittant = adder(emissions, emittants)

            pollutants_lst.append(pollutant)
            
            emittants_lst.append(sum_emittant)
            
        for i, month in enumerate(months):
            ax[i].add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
            ax[i].coastlines(resolution='50m', linewidth=0.5, color='black')
            
            print(len(emittants_lst[i]), len(pollutants_lst[i]))
            plt.scatter(emittants_lst, pollutants_lst)
            ax[i].set_title(f"{var}, monthly average {month}", fontsize=14)

        plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'division_processing', 'division_figures', f'{type_}_{var}_timeaveraged_Europe.png'))

end = time.time()
plt.show()

print(f"Process run in {end - start} s")
