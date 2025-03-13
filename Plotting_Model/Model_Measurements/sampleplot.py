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

# Import modules
import matplotlib.pyplot as plt  # For plotting
import cartopy.crs as ccrs  # For map projections
import cartopy.feature as cfeature  # For adding map features (coastlines, borders, etc.)
import xarray as xr  # For handling netCDF data


# Open DataSet and print an overview of it
ds = xr.open_dataset("PM25.annual.ON.nc4", engine="netcdf4")  # or use "h5netcdf"
print(ds)

# Select a DataArray
var = 'PM25'
da = ds[var]

# Select data from ground level at a specific time
# Select a region using 'lat' and the nearest 'lon' using the method='nearest'
        # Check the underlying numpy array

daSurf = da  # Keep the full latitude and longitude dimension
 # Nearest longitude
 # No 'lev' dimension to select

print(daSurf)                # Check the contents of the DataArray
print(daSurf.dtype)          # Check the data type of the DataArray
print(daSurf.values)         # Check the underlying numpy array



# Print the mean value of this selection
print(f'\nMean of {var} is = {float(daSurf.mean().values):.1f} {da.attrs["units"]}')

# Plot figure using matplotlib and cartopy. The use of some cosmetic options is illustrated here
fig = plt.figure(figsize=[12, 6])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')
daSurf.plot(transform=ccrs.PlateCarree(), vmin=0, vmax=50)
plt.show()