import xarray as xr
import os
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as ccrs

# Load dataset
ds = xr.open_dataset(os.path.join(os.path.dirname(__file__), "..", "raw_data", "misc", "lspop2019_05x0625.nc4"))
da = ds["lspop2019"].isel(lat = range(250, 320), lon = range(240, 368))


# Create figure and axis
fig, ax = plt.subplots(1, 1, figsize=[9, 7], subplot_kw={"projection": ccrs.EqualEarth(central_longitude=10)})
fig.tight_layout()

# Add map features
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='black')

# Plot data
da.plot(ax=ax,
        transform=ccrs.PlateCarree(),
        cmap= "hot_r"
        )

plt.show()
