# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import numpy as np
from utils import *
from sklearn.preprocessing import StandardScaler
import time
import seaborn as sns
from scipy.stats import zscore
from sklearn.svm import OneClassSVM
from sklearn.cluster import DBSCAN
start = time.time()

months = ['JAN', 'JUL']

types = {
    'Aerosol': {
        'PM25': ['nvPM'],
        'AerMassNIT': ['NO2'],
        'AerMassNH4': ['NO2'],
        'AerMassPOA': ['nvPM', 'HC'],
        'AerMassBC': ['nvPM']
    },
    'O3': {
        'SpeciesConc_O3': ['NO2', 'HC', 'CO']
    },
    'NO2': {
        'SpeciesConc_NO2': ['NO2']
    }
}
def modified_z_score(data):
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    return 0.6745 * (data - median) / mad


emissions = xr.open_dataset(os.path.join(os.path.dirname(__file__), '..', 'raw_data', 'emissions', 'AvEmMasses.nc4'))

# emittants:
# NO2, HC, CO, nvPM
# Iterate through each pollutant
for type_, vars in types.items():
    # Iterate through the emittants of each pollutant
    for var, emittants in vars.items():
        fig, axs = plt.subplots(1, 2, figsize=[12, 4])
        plt.tight_layout()
        measures = []
        
        for month in months:
            pollutant = differencer(type_, month, var)

            sum_emittants = adder(emissions, emittants)
            
            measure = pollutant/ (sum_emittants)     
            print(f'Minimum emittants {pollutant.min()} max emittants {pollutant.max()}')
            measures.append(measure)
            print(pollutant.coords)
            

        
        for i, month in enumerate(months):
            data = (measures[i].values.flatten())
            
            # # Calculate Q1 (25th percentile) and Q3 (75th percentile)
            # Q1 = np.percentile(data, 25)
            # Q3 = np.percentile(data, 75)

            # # Calculate the interquartile range (IQR)
            # IQR = Q3 - Q1
            
            # # Calculate the lower and upper bounds for outliers
            # lower_bound = Q1 - 25 * IQR
            # upper_bound = Q3 + 25 * IQR

            # # Identify outliers (values outside the whiskers)
            # outliers = data[(data < lower_bound) | (data > upper_bound)]

            # # Count outliers
            # num_outliers = len(outliers)
            
            # print(f"Number of outliers: {num_outliers}, data length {len(data)}")
            # sns.boxplot(data, ax=axs[i])
            axs[i].hist(data)
            
            
            
plt.show()   
            
end = time.time()


print(f"Process run in {end - start} s")