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

# ----------------------------------------------------
# 1. Prepare the land mask using Cartopy + Shapely
# ----------------------------------------------------
shpfilename = shpreader.natural_earth(resolution='10m', category='physical', name='land')
reader = shpreader.Reader(shpfilename)
geoms = reader.geometries()

polygons = []
for geom in geoms:
    if isinstance(geom, sgeom.MultiPolygon):
        polygons.extend(geom.geoms)
    else:
        polygons.append(geom)

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

# ----------------------------------------------------
# 3. Plot only over Europe
# ----------------------------------------------------
start = time.time()
months = ['JAN', 'JUL']

# Define bounds for Europe
lon_min, lon_max = -20, 30
lat_min, lat_max = 38, 60

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
        measures = []

        for month in months:
            pollutant = differencer(type_, month, var)
            sum_emittants = adder(emissions, emittants)

            # Slice pollutants, top 20% chosen
            pollutant_top20 = np.nanpercentile(pollutant.values, 70)
            select_pollutant = pollutant.where(pollutant.values >= pollutant_top20, np.nan)

            measure = select_pollutant  / sum_emittants

            
            measures.append(measure)

        for i, month in enumerate(months):
            ax[i].add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
            ax[i].coastlines(resolution='50m', linewidth=0.5, color='black')
            ax[i].set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

            masked_measure = mask_ocean(measures[i])
            masked_measure_europe = masked_measure.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
            
            
            # median = np.nanmedian(masked_measure_europe.values)
            # q1 = np.nanpercentile(masked_measure_europe.values, 25)
            # q3 = np.nanpercentile(masked_measure_europe.values, 75)
            # iqr = q3 - q1
            # lower_bound = q1 - 1.5 * iqr
            # upper_bound = q3 + 1.5 * iqr
            
            # # masked_measure_europe = masked_measure_europe.where((masked_measure_europe >= lower_bound) & (masked_measure_europe <= upper_bound), np.nan)
            
            vmin = np.nanmin(masked_measure_europe.values)
            vmax = np.nanmax(masked_measure_europe.values)
            
            
            lati, longi = np.unravel_index(np.nanargmax(masked_measure_europe.values), masked_measure_europe.shape)
            latitudes = masked_measure_europe['lat'].values
            longitudes = masked_measure_europe['lon'].values
            lat = latitudes[lati]
            lon = longitudes[longi]

            print(lat, lon)
            masked_measure_europe.plot(ax=ax[i], transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)
            ax[i].set_title(f"{var}, monthly average {month}", fontsize=14)

        plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'division_processing', 'division_figures', f'{type_}_{var}_timeaveraged_Europe.png'))

end = time.time()
plt.show()

print(f"Process run in {end - start} s")
