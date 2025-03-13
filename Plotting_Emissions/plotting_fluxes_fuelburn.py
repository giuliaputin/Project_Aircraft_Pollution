# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os

# Open DataSet and print an overview of it
ds = xr.open_dataset(os.path.join('raw_data', 'emissions', 'AvEmFluxes.nc4'))
print(ds)

# Select a DataArray
var = 'FUELBURN'
da = ds[var]

# Select data from ground level at a specific time
flight_level = 31  # Level index
daSurf = da.isel(lev=flight_level)  # Remove time selection if 'time' is not in the dataset

# Print dataset info
print(daSurf)
print("Min:", daSurf.min().values, "Max:", daSurf.max().values)
print("Mean:", daSurf.mean().values)

# **Apply scaling factor if necessary**
# Convert from kg/m²/s to nanograms/m²/s for better readability
scale_factor = 1e9  # Convert kg to nanograms (ng)
daSurf_scaled = daSurf * scale_factor

# Print mean in better format
mean_value = float(daSurf_scaled.mean().values)
print(f'\nMean of {var} is = {mean_value:.3f} ng/m²/s')

# **Adjust plot scale dynamically**
vmin = float(daSurf_scaled.min().values)
vmax = float(daSurf_scaled.max().values)

# Prevent vmax from being too small (adjust as needed)
if vmax < vmin + 1e-3:
    vmax = vmin + 1e-3  # Ensure there's enough contrast

# Plot figure using matplotlib and cartopy
fig = plt.figure(figsize=[12, 6])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Use scaled data to ensure visibility
daSurf_scaled.plot(transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)

# Update title to include flight level
plt.title(f'{var} (scaled to ng/m²/s) - Flight Level {flight_level}')
plt.show()

