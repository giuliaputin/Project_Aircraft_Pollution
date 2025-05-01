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

standard_density = 44.6 #mol/m^3
molar_mass_NO2 = 46.01 #g/mol
molar_mass_O3 = 48 #g/mol

#Define scaled index function
def assign_index_scaled(concentration):
    """
    Scales concentration to a 0–100 range based on percentiles,
    improving contrast while avoiding outlier distortion.
    """
    c = concentration.copy()
    min_val = np.nanpercentile(c, 1)
    max_val = np.nanpercentile(c, 99)
    scaled = (c - min_val) / (max_val - min_val)
    scaled = np.clip(scaled, 0, 1)
    return scaled * 100

# # Calculate index for each pollutant
# indexes = {}
# for pol in datasets:
#     indexes[pol] = assign_index_scaled(datasets[pol])
#
# # Combine indexes and scale again to highlight contrasts
# total_index = indexes['O3']*0.5 + indexes['NO2']*(25/60) + indexes['PM25']*(25/15)
# total_index = assign_index_scaled(total_index)

total_index = 0.5 * o3_diff.da_diff * standard_density * molar_mass_O3 * (1000000) + (25/60) *  no2_diff.da_diff * standard_density * molar_mass_NO2 * (1000000) + (25/15) * pm25_diff.da_diff

# Extract lat/lon
lats = datasets['O3']['lat'].values
lons = datasets['O3']['lon'].values
nlat, nlon = total_index.shape

# Find top 10 most polluted grid points
flat = total_index.values.flatten()
top10_idx = np.argpartition(flat, -1)[-1:]
top10_idx = top10_idx[np.argsort(flat[top10_idx])[::-1]]

# Get lat/lon and index values of those points
top10_coords = [(lats[i // nlon], lons[i % nlon], flat[i]) for i in top10_idx]

# Set global font size
plt.rcParams.update({'font.size': 40})  # You can increase this value for even larger text

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

# Add colorbar and title with larger font sizes
cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05)
cbar.set_label('Total Pollution Index (Aviation Contribution)', fontsize=23)
cbar.ax.tick_params(labelsize=20)

ax.set_title('Overall Aircraft Pollution Impact Over Europe in July', fontsize=25)

# Save figure
plt.savefig(r"C:\Users\unigi\Downloads\aircraft_pollution_europe.svg", format="svg", bbox_inches='tight')
plt.savefig("aircraft_pollution_europe_jul.pdf", format="pdf", bbox_inches='tight')

# Blue markers and labels
# for lat, lon, index_value in top10_coords:
#     ax.plot(lon, lat, 'bo', markersize=5, transform=ccrs.PlateCarree())
#     ax.text(lon + 0.5, lat + 0.5, f'{index_value:.1f}', color='blue', fontsize=10,
#             transform=ccrs.PlateCarree())
#     print(lat, lon, index_value)

plt.show()
