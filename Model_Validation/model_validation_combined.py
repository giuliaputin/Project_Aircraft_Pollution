import pandas as pd
import matplotlib.pyplot as plt
import os
import xarray as xr
import numpy as np

def process_measurements(file_path_csv, file_path_nc, model_var):
    # Load measurements
    df = pd.read_csv(file_path_csv)
    measurements = df[['Longitude', 'Latitude', 'Air Pollution Level']]

    # Load model data
    ds = xr.open_dataset(file_path_nc)
    lat = ds['lat'].values
    lon = ds['lon'].values
    model_data = ds[model_var].values

    model = pd.DataFrame({
        "Latitude": np.repeat(lat, len(lon)),
        "Longitude": np.tile(lon, len(lat)),
        "Air Pollution Level": model_data.flatten()
    }).dropna()

    # Filter measurement data based on quantiles
    lat_min, lat_max = measurements['Latitude'].quantile([0.01, 0.99])
    lon_min, lon_max = measurements['Longitude'].quantile([0.01, 0.99])
    measurements = measurements[
        (measurements['Latitude'] >= lat_min) & (measurements['Latitude'] <= lat_max) &
        (measurements['Longitude'] >= lon_min) & (measurements['Longitude'] <= lon_max)
    ]

    # Calculate difference
    validation = []
    for idx, row in measurements.iterrows():
        distance = ((model['Latitude'] - row['Latitude'])**2 + (model["Longitude"] - row["Longitude"])**2)**0.5
        closest_idx = distance.idxmin()
        validation.append(model['Air Pollution Level'].iloc[closest_idx] - row['Air Pollution Level'])

    measurements['difference'] = validation

    # Bias and bias-corrected normalized error
    bias = measurements['difference'].mean()
    NMB = measurements['difference'].sum() / measurements['Air Pollution Level'].sum()
    measurements['error'] = (measurements['difference'] - bias).abs() / measurements['Air Pollution Level'] * 100

    per_bias = bias / measurements['Air Pollution Level'].mean() * 100
    print(f"--- {model_var} ---")
    print(f'Percentage bias: {per_bias:.2f}%')
    print(f'NMB: {NMB:.4f}')
    print(f'Bias corrected normalized error: {measurements["error"].mean():.2f}%')

    # Remove extreme outliers (as in your NO2 script)
    threshold = 80000
    measurements = measurements[measurements['error'].abs() < threshold]

    return measurements['error']

# --- Process each pollutant ---
pm25_errors = process_measurements(
    file_path_csv="raw_data/measurements/PM25.csv",
    file_path_nc="Plotting_Model/Annual_Model_Measurements/PM25.annual.ON.nc4",
    model_var="PM25"
)

o3_errors = process_measurements(
    file_path_csv="raw_data/measurements/O3.csv",
    file_path_nc="raw_data/model/O3.annual.ON.nc4",
    model_var="SpeciesConc_O3"
)

no2_errors = process_measurements(
    file_path_csv="raw_data/measurements/NO2.csv",
    file_path_nc="raw_data/model/NO2.annual.ON.nc4",
    model_var="SpeciesConc_NO2"
)

# --- Plot all boxplots together ---
plt.figure(figsize=(8, 6))
plt.boxplot([pm25_errors, o3_errors, no2_errors], vert=True, showfliers=False, whis=[0, 100])
plt.xticks([1, 2, 3], ['PM2.5', 'O3', 'NO2'])
plt.yscale('log')
plt.ylabel('Bias-Corrected Normalized Error (%)')
plt.title('Bias-Corrected Normalized Error for PM2.5, O3, and NO2')
plt.grid(True)
plt.show()
