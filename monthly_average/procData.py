# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import numpy as np
from utils import differencer

import time
start = time.time()
months = ['JAN']

types = {
    'Aerosol': {
        'PM25': ['gino', 'giuseppe'],
        'AerMassNIT': [2],
        'AerMassNH4': [3],
        'AerMassPOA': [4],
        'AerMassBC': [5]
    },
    'O3': {
        'SpeciesConc_O3': [6]
    },
    'NO2': {
        'SpeciesConc_NO2': [7]
    }
}

# emittants:
# NO2, HC, CO, NVPM

    # Iterate through each pollutant
for type, vars in types.items():
    # Iterate through the emittants of each pollutant
    for var, emittants in vars.items():
        for emittant in emittants:
            for month in months:
                
                pollutant = differencer(type, month, var)
                
                
                
                
                
                




end = time.time()


print(f"Process run in {end - start} s")