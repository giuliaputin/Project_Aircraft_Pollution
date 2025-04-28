# Import necessary modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os

# Open the datasets
ds_on = xr.open_dataset(os.path.join(os.path.dirname(__file__), ".", 'NO2.JUL.ON.nc4'))   # Aviation ON
ds_off = xr.open_dataset(os.path.join(os.path.dirname(__file__), ".", 'NO2.JUL.ON.nc4'))  # Aviation OFF

# Select the variable for ozone concentration
var = 'SpeciesConc_NO2'  # Ozone variable
da_on = ds_on[var]
da_off = ds_off[var]

# Average over the vertical dimension (lev)
da_on_avg = da_on.mean(dim='lev')
da_off_avg = da_off.mean(dim='lev')

# Select the same time point for both datasets
time_point = '2019-01-15'
da_on_time = da_on_avg.sel(time=time_point, method='nearest')
da_off_time = da_off_avg.sel(time=time_point, method='nearest')

# Compute the difference (aviation contribution)
da_diff = da_on_time - da_off_time

# Manually set the units for da_diff
da_diff.attrs["units"] = da_on_time.attrs.get("units", "unknown")  # Use "unknown" if not found

# Print mean impact of aviation
print(f'\nMean aviation contribution to NO2 over Europe = {float(da_diff.mean().values):.1f} {da_diff.attrs["units"]}')

# Print mean impact of aviation


# Plot the difference
fig = plt.figure(figsize=[12, 8])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))  # European-centered projection

# Add features (borders, coastlines)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Plot the O3 difference using a diverging colormap
da_diff.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='RdBu_r', center=0,
             vmin=-abs(da_diff).max(), vmax=abs(da_diff).max())

# Set title and labels
ax.set_title('Aviation Contribution to Ozone (O₃) over Europe on 2019-01-15')

# Show the plot
plt.show()
