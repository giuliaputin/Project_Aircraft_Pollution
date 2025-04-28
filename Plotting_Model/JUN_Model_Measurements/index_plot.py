# Import necessary modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import numpy as np

# Pollutant settings
pollutants = {
    'O3': 'SpeciesConc_O3',
    'NO2': 'SpeciesConc_NO2',
    'PM25': 'SpeciesConc_PM25'
}

# Open the datasets
datasets = {}
for pol in pollutants.keys():
    datasets[pol+'_on'] = xr.open_dataset(f'{pol}.JUL.ON.nc4')
    datasets[pol+'_off'] = xr.open_dataset(f'{pol}.JUL.OFF.nc4')

# Choose the time point
time_point = '2019-01-15'

# Function to compute the aviation contribution
def compute_diff(ds_on, ds_off, varname):
    da_on = ds_on[varname]
    da_off = ds_off[varname]
    da_on_avg = da_on.mean(dim='lev')
    da_off_avg = da_off.mean(dim='lev')
    da_on_time = da_on_avg.sel(time=time_point, method='nearest')
    da_off_time = da_off_avg.sel(time=time_point, method='nearest')
    da_diff = da_on_time - da_off_time
    return da_diff

# Compute differences
diffs = {}
for pol, varname in pollutants.items():
    diffs[pol] = compute_diff(datasets[pol+'_on'], datasets[pol+'_off'], varname)

# Define function to assign index based on concentration
def assign_index(pollutant, concentration):
    if pollutant == 'NO2':
        bins = [0, 50, 100, 200, 400]
    elif pollutant == 'O3':
        bins = [0, 60, 120, 180, 240]
    elif pollutant == 'PM25':
        bins = [0, 15, 30, 55, 110]
    else:
        raise ValueError("Unknown pollutant")

    index = np.zeros_like(concentration)

    index = np.where(concentration <= bins[1], 12.5, index)          # 0-25
    index = np.where((concentration > bins[1]) & (concentration <= bins[2]), 37.5, index)   # 25-50
    index = np.where((concentration > bins[2]) & (concentration <= bins[3]), 62.5, index)   # 50-75
    index = np.where((concentration > bins[3]) & (concentration <= bins[4]), 87.5, index)   # 75-100
    index = np.where(concentration > bins[4], 125, index)             # >100

    return index

# Assign index for each pollutant
indexes = {}
for pol in ['O3', 'NO2', 'PM25']:
    indexes[pol] = assign_index(pol, diffs[pol])

# Sum the indexes
total_index = indexes['O3'] + indexes['NO2'] + indexes['PM25']

# Plot the result
fig = plt.figure(figsize=[12, 8])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))

# Add features (borders, coastlines)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Plot the total index
im = ax.pcolormesh(
    total_index['lon'], total_index['lat'], total_index,
    transform=ccrs.PlateCarree(),
    cmap='RdYlGn_r'
)

# Colorbar
cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05)
cbar.set_label('Total Pollution Index (Aviation Contribution)')

# Set title
ax.set_title('Overall Aviation Pollution Index over Europe on 2019-01-15')

# Show plot
plt.show()
