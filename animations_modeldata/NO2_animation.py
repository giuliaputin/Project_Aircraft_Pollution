# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import numpy as np
from difference_machine import differencer

# Select either July (JUL) or January (JAN)
month = 'JAN'

# select aviation ON or OFF
aviation = 'OFF'

# Select variable you want to animate
# Choose: ['SpeciesConc_NO2']
var = 'SpeciesConc_NO2'

# Select the level you want to animate [0, 72]
level = 1

# If average is set to True then the level variable above is overridden
# Average with respect to height is taken
average = True

subtract = True

# --------------------------------------------------------------------------------------------------------------------
# Starting of the preprocessing, no need to modify anything after this
# Open DataSet and print an overview of it
ds = xr.open_dataset(os.path.join(os.path.dirname(__file__), "..", 'raw_data', 'model', f'NO2.{month}.{aviation}.nc4'))

da = ds[var]


# Set up the figure and axis
fig = plt.figure(figsize=[12, 6])
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

for day in da.time.values:
    print(day)




# Function to update the plot for each time step
def update(frame):

    # Select the data for the current time step
    if subtract:
        daSurf = differencer("NO2",month,var,frame,level,average)
    else:
        if average:
            daSurf = da.mean(dim= "lev").isel(time=frame)
        else:
            daSurf = da.isel(lev = level).isel(time=frame)

    # Clear the previous plot
    ax.clear()
    
    # Redraw the coastlines and borders for every frame
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
    ax.coastlines(resolution='50m', linewidth=0.5, color='white')
    
    # Plot the data for the current time step
    im = daSurf.plot(ax=ax, transform=ccrs.PlateCarree(), add_colorbar=True)

    date_str = np.datetime_as_string(da.time[frame].values, unit='D') 

    if average:
        # Add a title with the current time
        ax.set_title(f"{var} at time {date_str}, averaged over height", fontsize=14)

    else:
        # Add a title with the current time
        ax.set_title(f"{var} at time {date_str}, level = {np.round(daSurf.lev.values, 3)}", fontsize=14)

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


