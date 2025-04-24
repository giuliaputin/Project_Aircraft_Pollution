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
from cartopy.io.shapereader import natural_earth
import shapely.geometry as sgeom
import cartopy.io.shapereader as shpreader
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.prepared import prep
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.prepared import prep

# ----------------------------------------------------
# 1. Prepare the land mask using Cartopy + Shapely
# ----------------------------------------------------

# Load natural earth land shapefile
shpfilename = shpreader.natural_earth(resolution='10m', category='physical', name='land')
reader = shpreader.Reader(shpfilename)
geoms = reader.geometries()

# Flatten MultiPolygons into list of Polygons
polygons = []
for geom in geoms:
    if isinstance(geom, sgeom.MultiPolygon):
        polygons.extend(geom.geoms)
    else:
        polygons.append(geom)

# Create a MultiPolygon and prepare it for fast masking
land_geom = sgeom.MultiPolygon(polygons)
prep_land = prep(land_geom)

# ----------------------------------------------------
# 2. Function to mask ocean values from a DataArray
# ----------------------------------------------------

def mask_ocean(dataarray):
    lons, lats = np.meshgrid(dataarray['lon'], dataarray['lat'])
    points = np.vstack([lons.ravel(), lats.ravel()]).T
    mask = np.array([prep_land.contains(sgeom.Point(x, y)) for x, y in points])
    mask = mask.reshape(lats.shape)
    return dataarray.where(mask)



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
        fig, ax = plt.subplots(1, 2, figsize= [12, 4], subplot_kw={"projection": ccrs.EqualEarth(central_longitude=10)})
        # fig, ax = plt.subplots(1, 1)
        plt.tight_layout()
        
        vmin = float('inf')
        vmax = float('-inf')
        measures = []
        
        for month in months:
            pollutant = differencer(type_, month, var)

            sum_emittants = adder(emissions, emittants)
            
            # sum_emittants = adder(emissions, emittants).drop_sel(lat=68.5).drop_isel(lon=-1)

            # This code below was not used
            # scaler = StandardScaler().fit(pollutant.values)
            
            # sum_emittants_scaled = xr.DataArray(
            #                         scaler.transform(sum_emittants.values),
            #                         dims=sum_emittants.dims,
            #                         coords=sum_emittants.coords,
            #                         attrs=sum_emittants.attrs
            #                     )
            
            # pollutant_scaled = xr.DataArray(
            #     scaler.transform(pollutant.values),
            #     dims=pollutant.dims,
            #     coords=pollutant.coords,
            #     attrs=pollutant.attrs
            # )
            #print(pollutant)
            measure = (pollutant)  / sum_emittants    
            
            median = np.median(measure.values)        
            q1 = np.percentile(measure.values, 25)
            q3 = np.percentile(measure.values, 75)
            iqr = q3 - q1

            lower_bound = q1 - 1.5 * iqr
            upper_bound = q3 + 1.5 * iqr

            # Old code
            # measure = measure.where(
            #     (measure >= lower_bound) & (measure <= upper_bound),
            #     median
            # )
            
            measure = measure.where(
                (measure >= lower_bound),
                np.nan
            )   
            measure = measure.where(
                (measure <= upper_bound),
                np.nan
            )
                 
            measures.append(np.log1p(np.abs(measure)))

        
        for i, month in enumerate(months):
            ax[i].add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
            ax[i].coastlines(resolution='50m', linewidth=0.5, color='black')

            # Mask out ocean values
            masked_measure = mask_ocean(measures[i])


            vmin = np.nanmin(masked_measure.values)
            vmax = np.nanmax(masked_measure.values)

            masked_measure.plot(ax=ax[i], transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)

            ax[i].set_title(f"{var}, monthly average {month}", fontsize=14)
            
            plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'division_processing', 'division_figures',f'{type_}_{var}_timeaveraged.png'))
end = time.time()
plt.show()

print(f"Process run in {end - start} s")