# Import libraries
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Instead of using xr.open_dataset, you import the python files
import O3_plotting_JUL_diff as o3_diff
import NO2_plotting_JUL_diff as no2_diff
import PM25_plotting_JUL_diff as pm25_diff

# Assume each .py file defines a variable called da_diff
pollutants = {
    'O3': o3_diff.da_diff,
    'NO2': no2_diff.da_diff,
    'PM25': pm25_diff.da_diff
}


# Open the datasets
for pol, (filename, varname) in pollutants.items():
    pollutants[pol] = xr.open_dataset(filename)[varname]


# Define the function to assign pollution index based on concentration
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

    index = np.where(concentration <= bins[1], 12.5, index)
    index = np.where((concentration > bins[1]) & (concentration <= bins[2]), 37.5, index)
    index = np.where((concentration > bins[2]) & (concentration <= bins[3]), 62.5, index)
    index = np.where((concentration > bins[3]) & (concentration <= bins[4]), 87.5, index)
    index = np.where(concentration > bins[4], 125, index)

    return index


# Assign index for each pollutant
indexes = {}
for pol in pollutants.keys():
    indexes[pol] = assign_index(pol, pollutants[pol])

# Sum the indexes
total_index = indexes['O3'] + indexes['NO2'] + indexes['PM25']

# Plotting
fig = plt.figure(figsize=(12, 8))
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))

# Add borders and coastlines
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Plot the total index
im = ax.pcolormesh(
    pollutants['O3']['lon'], pollutants['O3']['lat'], total_index,
    transform=ccrs.PlateCarree(),
    cmap='RdYlGn_r'
)

# Colorbar
cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05)
cbar.set_label('Total Pollution Index (Aviation Contribution)')

# Title
ax.set_title('Overall Aircraft Pollution Impact (July)')

# Show plot
plt.show()
