# Import libraries
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Import precomputed difference arrays
import O3_plotting_JUL_diff as o3_diff
import NO2_plotting_JUL_diff as no2_diff
import PM25_plotting_JUL_diff as pm25_diff

# Create dictionary of datasets
datasets = {
    'O3': o3_diff.da_diff,
    'NO2': no2_diff.da_diff,
    'PM25': pm25_diff.da_diff
}

# Define scaled index function
def assign_index_scaled(concentration):
    """
    Scales concentration to a 0–100 range based on percentiles,
    improving contrast while avoiding outlier distortion.
    """
    c = concentration.copy()
    min_val = np.nanpercentile(c, 5)
    max_val = np.nanpercentile(c, 95)
    scaled = (c - min_val) / (max_val - min_val)
    scaled = np.clip(scaled, 0, 1)
    return scaled * 100

# Calculate index for each pollutant
indexes = {}
for pol in datasets:
    indexes[pol] = assign_index_scaled(datasets[pol])

# Combine indexes and scale again to highlight contrasts
total_index = indexes['O3'] + indexes['NO2'] + indexes['PM25']
total_index = assign_index_scaled(total_index)

# Plotting
fig = plt.figure(figsize=(12, 8))
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Plot total pollution index
im = ax.pcolormesh(
    datasets['O3']['lon'], datasets['O3']['lat'], total_index,
    transform=ccrs.PlateCarree(),
    cmap='RdYlGn_r'
)

# Add colorbar and title
cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05)
cbar.set_label('Total Pollution Index (Aviation Contribution)')
ax.set_title('Overall Aircraft Pollution Impact (July)')

plt.show()
