import matplotlib.pyplot as plt
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import numpy as np
from utils import *
import time
from cartopy.io.shapereader import natural_earth
import shapely.geometry as sgeom
import cartopy.io.shapereader as shpreader
from shapely.prepared import prep
from matplotlib.patches import Rectangle

# Load the 'admin_0_countries' shapefile (countries of the world)
shpfilename = shpreader.natural_earth(resolution='10m', category='cultural', name='admin_0_countries')
reader = shpreader.Reader(shpfilename)
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
        if isinstance(geom, sgeom.MultiPolygon):
            polygons.extend(geom.geoms)
        else:
            polygons.append(geom)

# Create a MultiPolygon for Europe
europe_geom = sgeom.MultiPolygon(polygons)
prep_europe = prep(europe_geom)

# Rough bounding box for Europe
EUROPE_BOUNDS = {
    'lon_min': -25,
    'lon_max': 45,
    'lat_min': 30,
    'lat_max': 75
}

# ----------------------------------------------------
# Function to mask ocean values from a DataArray
# ----------------------------------------------------
def mask_ocean(dataarray):
    lons, lats = np.meshgrid(dataarray['lon'], dataarray['lat'])
    points = np.vstack([lons.ravel(), lats.ravel()]).T
    mask = np.array([prep_europe.contains(sgeom.Point(x, y)) for x, y in points])
    mask = mask.reshape(lats.shape)
    return dataarray.where(mask)

# ----------------------------------------------------
# Plot only over Europe
# ----------------------------------------------------
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
        
        measures = []
        
        for month in months:
            pollutant = differencer(type_, month, var)
            sum_emittants = adder(emissions, emittants)

            pollutant_top20 = np.nanpercentile(pollutant.values, 100*0)
            select_pollutant = pollutant.where(pollutant.values >= pollutant_top20, np.nan)

            measure = select_pollutant / sum_emittants
            
            measures.append(measure)
            # 48.125 in lon cutoff
            # 68.5 in lat is cutoff
        
        fig, ax = plt.subplots(1, 2, figsize=[16, 6], subplot_kw={"projection": ccrs.PlateCarree()})
        plt.tight_layout()
        
        for i, month in enumerate(months):
            ax[i].add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=0.5, edgecolor='darkgrey')
            ax[i].coastlines(resolution='10m', linewidth=0.5, color='black')
            
            # Set extent to Europe
            ax[i].set_extent([EUROPE_BOUNDS['lon_min'], EUROPE_BOUNDS['lon_max'],
                              EUROPE_BOUNDS['lat_min'], EUROPE_BOUNDS['lat_max']], crs=ccrs.PlateCarree())
            
            # Add bounding box on map
            europe_rect = Rectangle(
                (EUROPE_BOUNDS['lon_min'], EUROPE_BOUNDS['lat_min']),
                EUROPE_BOUNDS['lon_max'] - EUROPE_BOUNDS['lon_min'],
                EUROPE_BOUNDS['lat_max'] - EUROPE_BOUNDS['lat_min'],
                linewidth=1.5, edgecolor='blue', facecolor='none', linestyle='--', transform=ccrs.PlateCarree()
            )
            ax[i].add_patch(europe_rect)

            # Optional: Add gridlines with labels
            gl = ax[i].gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5, linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {'size': 9}
            gl.ylabel_style = {'size': 9}

            masked_measure = mask_ocean(measures[i])
            negative_measures = abs(masked_measure.where(masked_measure.values < 0))
            masked_measure_europe = masked_measure.where(masked_measure.values > 0)
            
            if IQR:
                median = np.nanmedian(masked_measure_europe.values)
                q1 = np.nanpercentile(masked_measure_europe.values, 25)
                q3 = np.nanpercentile(masked_measure_europe.values, 75)
                iqr = q3 - q1
                lower_bound = q1 - 1.5 * iqr
                upper_bound = q3 + 1.5 * iqr
                masked_measure_europe = masked_measure_europe.where(
                    (masked_measure_europe >= lower_bound) & 
                    (masked_measure_europe <= upper_bound), np.nan
                )

            vmin = np.nanmin(masked_measure_europe.values)
            vmax = np.nanmax(masked_measure_europe.values)
            norm1 = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
            
            fig.set_size_inches(16, 9)
            if np.any(~np.isnan(negative_measures.values)):
                norm2 = matplotlib.colors.LogNorm(vmin=np.nanmin(negative_measures.values), vmax=np.nanmax(negative_measures.values))
                negative_measures.plot(
                    ax=ax[i],
                    cmap='BuGn',
                    transform=ccrs.PlateCarree(),
                    norm=norm2,
                    cbar_kwargs={"label": "Negative values (log scale)", "aspect": 40, "location": "bottom", "pad": 0.05}
                )
            
            masked_measure_europe.plot(
                ax=ax[i],
                transform=ccrs.PlateCarree(),
                norm=norm1,
                cmap='YlOrRd',
                cbar_kwargs={"label": "Absolute value (log scale)", "aspect": 40, "location": "bottom", "pad": 0.05}
            )
            
            ax[i].set_title(f"{var}, monthly average {month}", fontsize=14)

        plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'division_processing', 'division_figures', f'{type_}_{var}_timeaveraged_logscale.png'))

end = time.time()
# plt.show()

print(f"Process run in {end - start} s")
