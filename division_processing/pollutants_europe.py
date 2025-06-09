# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import os
from utils import differencer
# Types and variables that will be plotted
types = {
    "Aerosol": {
        "PM25": {
            "emittants": ["NO2"],
            "unit": 1,
            "name": r'PM$_{2.5} \:$'
        },
        "AerMassNIT": {
            "emittants": ["NO2"],
            "unit": 1,  # Example, adjust as needed
            "name": r'Inorganic Nitrate Aerosols'
        },
        "AerMassSO4":{
            "emittants": ["FUELBURN"],
            "unit": 1,
            "name": r"Sulfate Aerosol"
        }
        # 'AerMassNH4': ['NO2'],
        # 'AerMassPOA': ['nvPM', 'HC'],
        # 'AerMassBC': ['nvPM']
    },
    # "Aerosol2": {
    #     "AerMassNIT": ["NO2"]
    # },
    "O3": {
        "SpeciesConc_O3": {
            "emittants": ["NO2", "HC", "CO"],
            "unit": 44.6 * 48 * 1e6,
            "name": r'O$_3 \:$'
        }
    },
    "NO2": {
        "SpeciesConc_NO2": {
            "emittants": ["NO2"],
            "unit": 44.6 * 46.01 * 1e6,
            "name": r'NO$_2 \:$'
        }
    },
}

# --------------------------------------------------------------------------------------------------------------------
# Starting of the preprocessing, no need to modify anything after this

months = ['JAN', 'JUL']

# For loop to iterate through each type and through each variable withing each type
for type_, vars in types.items():
    for var, properties in vars.items():
        unit = properties["unit"]
        name = properties["name"]
        fig, ax = plt.subplots(1, 2, figsize= [12, 4], subplot_kw={"projection": ccrs.EqualEarth(central_longitude=10)})
        plt.tight_layout()

        daSurf = [[],[]]
        for i, month in enumerate(months):

            daSurf[i] = differencer(type_, month, var) * unit

            if i==0:
                vmin = float(daSurf[i].values.min())
                vmax = float(daSurf[i].values.max())
            else:
                vmin = min(vmin, float(daSurf[i].values.min()))
                vmax = max(vmax, float(daSurf[i].values.max()))

        for i, month in enumerate(months):

            ax[i].add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
            ax[i].coastlines(resolution='50m', linewidth=0.5, color='black')

            # Plot the data for the current time step
            daSurf[i].plot(ax=ax[i], transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cbar_kwargs ={"label": name +r" [$\mu g /  m^3$]"})

            # Add a title with the current time
            ax[i].set_title(month, fontsize=14)
        
       
        plt.savefig(os.path.join(os.path.dirname(__file__), '..', 'division_processing', 'pollutants_figures',f'{type_}_{var}_timeaveraged.png'))
        
# Uncomment to display the plots on screen
# plt.show()