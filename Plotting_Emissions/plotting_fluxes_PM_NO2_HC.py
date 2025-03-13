# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr

# Open DataSet and print an overview of it
ds = xr.open_dataset('raw_data/emissions/AvemFluxes.nc4')
print(ds)

# Select a DataArray
var = 'nvPM'  
da = ds[var]

# Select data from ground level
daSurf = da.isel(lev=30)  # Level 0 selection

# Print dataset info
print(daSurf)
print("Min:", daSurf.min().values, "Max:", daSurf.max().values)
print("Mean:", daSurf.mean().values)

# **Apply scaling factor for better readability**
scale_factor = 1e15  # Convert from kg/m²/s to picograms per m²/s (pg/m²/s)
daSurf_scaled = daSurf * scale_factor

# Compute the new mean after scaling
mean_value = float(daSurf_scaled.mean().values)
print(f'\nMean of {var} is = {mean_value:.3f} pg/m²/s')

# **Adjust plot scale dynamically**
vmin = float(daSurf_scaled.min().values)
vmax = float(daSurf_scaled.max().values)

# Prevent vmax from being too small (adjust as needed)
if vmax < vmin + 1e-3:
    vmax = vmin + 1e-3  # Ensure there's enough contrast

# **Plot the scaled data**
fig = plt.figure(figsize=[12, 6])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Use scaled data to ensure visibility
daSurf_scaled.plot(transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)

plt.title(f'{var} (scaled to pg/m²/s)')
plt.show()
