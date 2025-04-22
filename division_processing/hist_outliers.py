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
import seaborn as sns
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

# emittants:
# NO2, HC, CO, nvPM
# Iterate through each pollutant
for type_, vars in types.items():
    # Iterate through the emittants of each pollutant
    for var, emittants in vars.items():
        # fig, ax = plt.subplots(1, 1, figsize=[12, 4])
        plt.tight_layout()
        measures = []
        
        for month in months:
            pollutant = differencer(type_, month, var)

            sum_emittants = adder(emissions, emittants).drop_sel(lat=68.5).drop_isel(lon=-1)

            scaler = StandardScaler().fit(pollutant.values)
            
            sum_emittants_scaled = xr.DataArray(
                                    scaler.transform(sum_emittants.values),
                                    dims=sum_emittants.dims,
                                    coords=sum_emittants.coords,
                                    attrs=sum_emittants.attrs
                                )
            
            pollutant_scaled = xr.DataArray(
                scaler.transform(pollutant.values),
                dims=pollutant.dims,
                coords=pollutant.coords,
                attrs=pollutant.attrs
            )
            
            measure = (pollutant) / sum_emittants    
        
            # Old code
            # measure = measure.where(
            #     (measure >= lower_bound) & (measure <= upper_bound),
            #     median
            # )
        
            measures.append(measure)

        
        for i, month in enumerate(months):
        
            
            # Plot the data for the current time step
            sns.boxplot(x= measures[i].values.flatten())
            # Add a title with the current time
            # ax.set_title(f"{var}, monthly average {month}", fontsize=14)
            plt.show()
            # plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'division_processing', 'division_figures',f'{type_}_{var}_timeaveraged.png'))
end = time.time()


print(f"Process run in {end - start} s")