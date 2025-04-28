# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os

# Open the dataset and print an overview
ds = xr.open_dataset(os.path.join(os.path.dirname(__file__), "..", "..", 'raw_data', 'model', 'NO2.JUL.OFF.nc4'))
print(ds)

# Select the variable for ozone concentration
var = 'SpeciesConc_NO2'  # Corrected variable name for ozone
da = ds[var]

# Average over the vertical dimension (lev) to get the total or averaged O3 concentration at each location
da_avg = da.mean(dim='lev')

# Select data at a specific time, using nearest available time point
da_time = da_avg.sel(time='2019-06-15', method='nearest')

# Print the mean value of this selection
print(f'\nMean of {var} over Europe is = {float(da_time.mean().values):.1f} {da.attrs["units"]}')

# Plot the spatial distribution of O3 over Europe
fig = plt.figure(figsize=[12, 8])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))  # European-centered projection

# Add features to the map (borders, coastlines)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Plot the O3 data using a better colormap and adjust the color limits dynamically
da_time.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='viridis',
             vmin=da_time.min().values, vmax=da_time.max().values)

# Set title and labels
ax.set_title('Ozone (O₃) Distribution over Europe on 2019-01-15')

# Show the plot
plt.show()
