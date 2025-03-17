import matplotlib.pyplot as plt
import os
import sys
import numpy as np
import xarray as xr #netCDF Library
import cartopy
import cartopy.crs as ccrs #projections list
import cartopy.feature as cfeature

# for venv use
print(sys.executable)

def differencer(type, month, var, frame, level, average=False):
    dsoff = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." ,"raw_data", "model", f"{type}.{month}.OFF.nc4") )
    dson = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." ,"raw_data", "model", f"{type}.{month}.ON.nc4") )
    daoff = dsoff[var]
    daon = dson[var]
    if average:
        daSurfoff = daoff.isel(time=frame).mean(dim = "lev")
        daSurfon = daon.isel(time=frame).mean(dim = "lev")
    else:
        daSurfoff = daoff.isel(time=frame).isel(lev=level)
        daSurfon = daon.isel(time=frame).isel(lev=level)
    daDiff = daSurfon - daSurfoff

    return daDiff


# Some defaults:

plt.rcParams['figure.figsize'] = (12, 5)  # Default plot size

# dsoff = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." ,"raw_data", "model", "Aerosol.JAN.OFF.nc4") )
# dson = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." ,"raw_data", "model", "Aerosol.JAN.ON.nc4") )
    

# offPM = dsoff["PM25"]
# onPM = dson["PM25"]


# #print(offPM.isel(time=0))
# altsumoff = offPM.isel(time=0).isel(lev=71)
# altsumon = onPM.isel(time=0).isel(lev=71)

# print("Bananas")

# diff = altsumon - altsumoff

# diff = differencer("Aerosol", "JAN", "PM25", 0, 71)


# # Define the map projection (how does the map look like)
# ax = plt.axes(projection=ccrs.EqualEarth())
# # ax is an empty plot. We now plot the variable sw_avg onto ax
# diff.plot(ax=ax, transform=ccrs.PlateCarree()) 
# # the keyword "transform" tells the function in which projection the data is stored 
# ax.coastlines(); ax.gridlines(); ax.add_feature(cfeature.BORDERS.with_scale('50m')); # Add gridlines and coastlines to the plot
# plt.show()

# # Define the map projection (how does the map look like)
# ax = plt.axes(projection=ccrs.EqualEarth())
# # ax is an empty plot. We now plot the variable sw_avg onto ax
# altsum.plot(ax=ax, transform=ccrs.PlateCarree()) 
# # the keyword "transform" tells the function in which projection the data is stored 
# ax.coastlines(); ax.gridlines(); # Add gridlines and coastlines to the plot