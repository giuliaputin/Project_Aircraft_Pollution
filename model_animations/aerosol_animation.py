# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from difference_machine import differencer

# Select either July (JUL) or January (JAN)
month = 'JAN'

# select aviation ON or OFF
aviation = 'ON'

# Select variable you want to animate
# Choose between: ['PM25', 'AerMassNIT', 'AerMassNH4', 'AerMassPOA', 'AerMassBC']
var = 'PM25'

level = 1
# --------------------------------------------------------------------------------------------------------------------
# Starting of the preprocessing, no need to modify anything after this
# Open DataSet and print an overview of it
ds = xr.open_dataset(os.path.join(os.path.dirname(__file__),"..",'raw_data', 'model', f'Aerosol.{month}.{aviation}.nc4'))

da = ds[var]

# Set up the figure and axis
fig = plt.figure(figsize=[12, 6])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

average = False
# Function to update the plot for each time step
def update(frame):

    # Select the data for the current time step
    
    daSurf = differencer("Aerosol",month,var,frame,level,average)
    
    # Clear the previous plot
    ax.clear()
    
    # Redraw the coastlines and borders for every frame
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
    ax.coastlines(resolution='50m', linewidth=0.5, color='white')
    
    # Plot the data for the current time step
    im = daSurf.plot(ax=ax, transform=ccrs.PlateCarree(), add_colorbar=False)
    #vmin=0, vmax=50
    # Add a title with the current time
    ax.set_title(f"{var} at time {da.time[frame].values}", fontsize=14)
    
    return [im]

# Create the animation
ani = FuncAnimation(fig, update, frames=len(da.time), interval=500, blit=False)

# Display the animation
plt.show()

# ani.save(os.path.join('model_animations', 'animations', f'{month}_{aviation}_{var}.mp4'), writer='ffmpeg', fps=4)
