# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from difference_machine import differencer
import numpy as np

# From here you can set up the animation:
# Select either July (JUL) or January (JAN)
month = 'JAN'

# select aviation ON or OFF
aviation = 'OFF'

# Select variable you want to animate
# Choose between: ['PM25', 'AerMassNIT', 'AerMassNH4', 'AerMassPOA', 'AerMassBC']
var = 'AerMassBC'

# Select the level you want to animate [0, 72]
level = 0

# If average is set to True then the level variable above is overridden
# Average with respect to height is taken
average = True

# select subract to find aviation specific emissions. This overrides the aviation = "on"
subtract = False

# --------------------------------------------------------------------------------------------------------------------
# Starting of the preprocessing, no need to modify anything after this
# Open DataSet and print an overview of it
ds = xr.open_dataset(os.path.join(os.path.dirname(__file__),"..",'raw_data', 'model', f'Aerosol.{month}.{aviation}.nc4'))

da = ds[var]

# Set up the figure and axis
fig = plt.figure(figsize=[12, 6])
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

vmin, vmax = np.inf, -np.inf


for frame in range(len(da.time)):
    if subtract:
        daSurf = differencer("Aerosol", month, var, frame, level, average)
    else:
        if average:
            daSurf = da.mean(dim="lev").isel(time=frame)
        else:
            daSurf = da.isel(lev=level).isel(time=frame)
    vmin = min(vmin, float(daSurf.min(skipna=True)))
    vmax = max(vmax, float(daSurf.max(skipna=True)))

print(f"Global vmin: {vmin}, vmax: {vmax}")

# print(vmin, vmax)
# Function to update the plot for each time step
def update(frame):
    # Select the data for the current time step
    if subtract:
        daSurf = differencer("Aerosol", month, var, frame, level, average)
    else:
        if average:
            daSurf = da.mean(dim="lev").isel(time=frame)
        else:
            daSurf = da.isel(lev=level).isel(time=frame)
    
    # Normalize the data between 0 and 1
    min_val = float(daSurf.min(skipna=True))
    max_val = float(daSurf.max(skipna=True))
    daSurf_norm = (daSurf - min_val) / (max_val - min_val + 1e-10)  # Add epsilon to avoid division by zero
    
    # Clear the axis to prevent overplotting
    ax.clear()
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
    ax.coastlines(resolution='50m', linewidth=0.5, color='white')
    
    # Plot the normalized data
    im = daSurf_norm.plot(ax=ax, transform=ccrs.PlateCarree(), add_colorbar=False, cmap="viridis")
    
    date_str = np.datetime_as_string(da.time[frame].values, unit='D')
    if average:
        ax.set_title(f"{var} (normalized) at time {date_str}, averaged over height", fontsize=14)
    else:
        ax.set_title(f"{var} (normalized) at time {date_str}, level = {np.round(daSurf.lev.values, 3)}", fontsize=14)

    return [im]


# Create the animation
ani = FuncAnimation(fig, update, frames=len(da.time), interval=500, blit=False)

# Display the animation

plt.show()

# Saving of the files
if average:
    ani.save(os.path.join(os.path.dirname(__file__),"..",'model_animations', 'animations', f'{month}_{aviation}_{var}_avg.mp4'), writer='ffmpeg', fps=4)
else:
    ani.save(os.path.join(os.path.dirname(__file__),"..", 'model_animations', 'animations', f'{month}_{aviation}_{var}_lev{level}.mp4'), writer='ffmpeg', fps=4)

