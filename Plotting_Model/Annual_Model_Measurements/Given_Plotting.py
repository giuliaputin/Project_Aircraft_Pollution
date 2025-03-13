# sample_plot.py
#
# This script is an example showing how you can use [xarray] to access the
#  data provided in the netCDF files and then use [matplotlib] and [cartopy]
#  to plot a figure using the data
#
# For more information on these modules see:
#  https://docs.xarray.dev/en/stable/
#  https://matplotlib.org/index.html
#  https://scitools.org.uk/cartopy/docs/latest/
#
# — Flávio Quadros, 05/02/2020
#    - updated 07/03/2024

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr

# Open DataSet and print an overview of it
ds = xr.open_dataset('PM25.annual.ON.nc4')
print(ds)

# Select a DataArray
var = 'PM25'
da = ds[var]

# Since there's no "lev" or "time" dimension, use the dataset directly
daSurf = da

# Print the mean value of this selection
print(f'\nMean of {var} is = {float(daSurf.mean().values):.1f}')

# Plot figure using matplotlib and cartopy
fig = plt.figure(figsize=[12, 6])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Plot the data
daSurf.plot(transform=ccrs.PlateCarree(), vmin=0, vmax=50)
plt.show()
