import matplotlib.pyplot as plt
import xarray as xr
from shapely.geometry import MultiPolygon
from shapely.vectorized import contains
import os
import numpy as np
from utils import *
import time
from cartopy.io.shapereader import natural_earth, Reader
from shapely.prepared import prep

start = time.time()
months = ['JAN', 'JUL']


# Load the 'admin_0_countries' shapefile (countries of the world)
shpfilename = natural_earth(resolution='10m', category='cultural', name='admin_0_countries')
reader = Reader(shpfilename)
records = reader.records()

IQR = False
europe_countries = [
    "Austria", "Belgium", "Bosnia and Herzegovina", "Bulgaria", "Croatia", "Czech Republic",
    "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary",
    "Iceland", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta",
    "Moldova", "Netherlands", "North Macedonia", "Norway", "Poland", "Portugal",
    "Romania", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden", "Switzerland",
    "Ukraine", "United Kingdom", "Albania", "Montenegro", "Kosovo", "Belarus"
]

# Collect geometries for only the European countries
polygons = []
for record in records:
    if record.attributes['NAME_LONG'] in europe_countries:
        geom = record.geometry
        if isinstance(geom, MultiPolygon):
            polygons.extend(geom.geoms)
        else:
            polygons.append(geom)

# Create a MultiPolygon for Europe
europe_geom = MultiPolygon(polygons)
prep_europe = prep(europe_geom)

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
        fig, ax = plt.subplots(1, 2, figsize=[12, 3])
        plt.tight_layout()
        pollutants_lst = []
        emittants_lst = []
        
        for month in months:
            pollutant = differencer(type_, month, var) # 74 by 122
            sum_emittant = adder(emissions, emittants) # 75 by 123
            
            sum_emittant = (sum_emittant.drop_isel({'lat':-1, 'lon': -1}))
            
            # Filter data for europe 
            lon2d, lat2d = np.meshgrid(sum_emittant['lon'].values, sum_emittant['lat'].values)
            mask = contains(europe_geom, lon2d, lat2d)
            
            pollutants_lst.append(pollutant.where(mask))
            
            emittants_lst.append(sum_emittant.where(mask)) 
            
            
        for i, month in enumerate(months):
            np.set_printoptions(threshold=np.inf)

            # mask the top emissions data
            threshold = np.nanpercentile(emittants_lst[i], 97.5)
            mask_emit = emittants_lst[i] > threshold
            
            ax[i].scatter(emittants_lst[i].where(~mask_emit), pollutants_lst[i].where(~mask_emit), color='blue')
            ax[i].scatter(emittants_lst[i].where(mask_emit), pollutants_lst[i].where(mask_emit), color='orange')
            ax[i].set_title(f"{var}, monthly average {month}", fontsize=14)
            ax[i].set_xscale('log')
            ax[i].set_yscale('log')
            ax[i].set_xlabel('Sum Emittants (kg/year)')
            ax[i].set_ylabel('Pollutants')
        plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'division_processing', 'scatter_figures_log', f'{type_}_{var}_timeaveraged_Europe_scatter_log.png'))

end = time.time()
plt.show()

print(f"Process run in {end - start} s")
