# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import numpy as np


# Open DataSet and print an overview of it
ds_emission = xr.open_dataset(os.path.join(os.path.dirname(__file__),"..",'raw_data', 'emissions', 'AvEmMasses.nc4'))

# Select a DataArray
vars = ["FUELBURN", "NO2", "HC", "CO", 'nvPM']
da_emission = ds_emission[vars[1]]
# for var in vars[1:]:
#     da_add = ds_emission[var]
#     da_emission += da_add

print(da_emission)

# fig = plt.figure(figsize=[12, 4])
# plt.tight_layout()
# ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
# ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
# ax.coastlines(resolution='50m', linewidth=0.5, color='black')
# ds_emission.plot(ax=ax, transform=ccrs.PlateCarree())

# plt.show()

# # Types and variables that will be plotted
# types = {'Aerosol': ['PM25', 'AerMassNIT', 'AerMassNH4', 'AerMassPOA', 'AerMassBC'],
#                   'O3' : ['SpeciesConc_O3'],
#                   'NO2' : ['SpeciesConc_NO2']}
               

# def differencer(type, month, var):
#     dsoff = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." ,"raw_data", "model", f"{type}.{month}.OFF.nc4") )
#     dson = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." ,"raw_data", "model", f"{type}.{month}.ON.nc4") )
#     daoff = dsoff[var]
#     daon = dson[var]

#     daSurfoff = daoff.isel(lev=0).mean(dim = "time")
#     daSurfon = daon.isel(lev=0).mean(dim = "time")
#     daDiff = daSurfon - daSurfoff

#     return daDiff

# # --------------------------------------------------------------------------------------------------------------------
# # Starting of the preprocessing, no need to modify anything after this

# months = ['Jan', 'JUL']

# # For loop to iterate through each type and through each variable withing each type
# for type in types:
#     for var in types[type]:
#         fig, ax = plt.subplots(1, 2, figsize= [12, 4], subplot_kw={"projection": ccrs.EqualEarth(central_longitude=10)})
#         plt.tight_layout()
    
#         for i, month in enumerate(months):

#             daSurf = differencer(type, month, var)
#             ax[i].add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
#             ax[i].coastlines(resolution='50m', linewidth=0.5, color='black')

#             # Plot the data for the current time step
#             daSurf.plot(ax=ax[i], transform=ccrs.PlateCarree(), add_colorbar=False)

#             # Add a title with the current time
#             ax[i].set_title(f"{var}, monthly average {month}, level = {np.round(daSurf.lev.values, 3)}", fontsize=14)
        
       
#         plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'monthly_average', 'figures',f'{type}_{var}_timeaveraged.png'))
        
# # Uncomment to display the plots on screen
# # plt.show()