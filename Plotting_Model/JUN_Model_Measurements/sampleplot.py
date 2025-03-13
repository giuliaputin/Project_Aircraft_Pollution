# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr

# Open DataSet and print an overview of it
ds = xr.open_dataset('O3.JUL.ON.nc4')
print(ds)

# Select a DataArray
var = 'SpeciesConc_O3'  # Corrected variable name for ozone
da = ds[var]

# Select data from ground level at the closest time to 2019-01-15
daSurf = da.isel(lev=0).sel(time='2019-01-15', method='nearest')

# Print the mean value of this selection
print(f'\nMean of {var} is = {float(daSurf.mean().values):.1f} {da.attrs["units"]}')

# Plot figure using matplotlib and cartopy. The use of some cosmetic options is illustrated here
fig = plt.figure(figsize=[12, 6])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')
daSurf.plot(transform=ccrs.PlateCarree(), vmin=0, vmax=50)
plt.show()
