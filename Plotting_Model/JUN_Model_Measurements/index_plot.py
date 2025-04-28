# Import necessary modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import numpy as np


# Function to assign index based on pollutant concentration
def assign_index(pollutant, concentration):
    if pollutant == 'NO2':
        if concentration <= 50:
            return 12.5
        elif concentration <= 100:
            return 37.5
        elif concentration <= 200:
            return 62.5
        elif concentration <= 400:
            return 87.5
        else:
            return 125
    elif pollutant == 'O3':
        if concentration <= 60:
            return 12.5
        elif concentration <= 120:
            return 37.5
        elif concentration <= 180:
            return 62.5
        elif concentration <= 240:
            return 87.5
        else:
            return 125
    elif pollutant == 'PM25':
        if concentration <= 15:
            return 12.5
        elif concentration <= 30:
            return 37.5
        elif concentration <= 55:
            return 62.5
        elif concentration <= 110:
            return 87.5
        else:
            return 125
    else:
        raise ValueError('Unknown pollutant')


# Open the datasets
datasets = {}
variables = {'NO2': 'SpeciesConc_NO2', 'O3': 'SpeciesConc_O3', 'PM25': 'SpeciesConc_PM25'}

fodatasets = {}

# Explicit filenames
datasets['NO2_on'] = xr.open_dataset('your_filename_for_NO2_ON.nc4')
datasets['NO2_off'] = xr.open_dataset('your_filename_for_NO2_OFF.nc4')

datasets['O3_on'] = xr.open_dataset('O3.JUL.ON.nc4')    # If this is correct
datasets['O3_off'] = xr.open_dataset('O3.JUL.OFF.nc4')

datasets['PM25_on'] = xr.open_dataset('your_filename_for_PM25_ON.nc4')
datasets['PM25_off'] = xr.open_dataset('your_filename_for_PM25_OFF.nc4')


# Select the time point
time_point = '2019-01-15'

# Calculate the aviation contribution and assign indexes
indexes = []

for pollutant in ['NO2', 'O3', 'PM25']:
    var = variables[pollutant]

    da_on = datasets[pollutant + '_on'][var].mean(dim='lev').sel(time=time_point, method='nearest')
    da_off = datasets[pollutant + '_off'][var].mean(dim='lev').sel(time=time_point, method='nearest')

    da_diff = da_on - da_off

    # Apply index assignment function vectorized
    index = xr.apply_ufunc(lambda x: assign_index(pollutant, x),
                           da_diff,
                           vectorize=True)

    indexes.append(index)

# Sum the three indexes
total_index = sum(indexes)

# Plot the total index
fig = plt.figure(figsize=[12, 8])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))

# Add features (borders, coastlines)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Plot total index
cmap = plt.get_cmap('Reds')
total_index.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap,
                 vmin=0, vmax=375)  # 3 pollutants * 125 max index

ax.set_title('Total Aviation Pollution Impact Index over Europe on 2019-01-15')
plt.show()
