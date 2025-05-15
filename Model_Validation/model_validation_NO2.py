import pandas as pd
import matplotlib.pyplot as plt
import folium
from folium.plugins import HeatMap
import os
import xarray as xr
import numpy as np

# Define the file path
file_path = "raw_data/measurements/NO2.csv"  # Adjust based on your folder structure

# Check if the file exists
if not os.path.exists(file_path):
    raise FileNotFoundError(f"File not found: {file_path}. Check your folder structure!")

# Load the CSV file
df = pd.read_csv(file_path)

# Check if required columns exist
required_columns = {"Longitude", "Latitude", "Air Pollution Level"}
if not required_columns.issubset(df.columns):
    raise ValueError(f"CSV file is missing one or more required columns: {required_columns}")




measurements = df[['Longitude', 'Latitude', 'Air Pollution Level']]

# Load the dataset
file_path = "raw_data/model/NO2.annual.ON.nc4"
ds = xr.open_dataset(file_path)

# Extract ongitude, and PM2.5 data
lat = ds['lat'].values
lon = ds['lon'].values
NO2 = ds['SpeciesConc_NO2'].values

# Convert to Pandas DataFrame
model = pd.DataFrame({
    "Latitude": np.repeat(lat, len(lon)),  # Repeat lat values
    "Longitude": np.tile(lon, len(lat)),   # Tile lon values
    "NO2": NO2.flatten()                 # Flatten NO2 values
})

# Remove NaN values
model = model.dropna()

model = model.rename(columns={'NO2':'Air Pollution Level'})
# Display the full data in table form (save to CSV if needed)
model.to_csv("no2_data.csv", index=False)

# model = model[(model['distance'] <= 1)]
# Define bounds (adjust if needed)
lat_min, lat_max = measurements['Latitude'].quantile([0.01, 0.99])
lon_min, lon_max = measurements['Longitude'].quantile([0.01, 0.99])

# Filter DataFrame
measurements = measurements[
    (measurements['Latitude'] >= lat_min) & (measurements['Latitude'] <= lat_max) &
    (measurements['Longitude'] >= lon_min) & (measurements['Longitude'] <= lon_max)
]



#uiteindelijke df

validation = []
distance = []
for idx, row in measurements.iterrows():
    validation.append(model['Air Pollution Level'].iloc[(((model['Latitude'] - row['Latitude'])**2 + (model["Longitude"] - row["Longitude"])**2)**0.5).idxmin()] - row['Air Pollution Level'])
    # distance.append((((model['Latitude'] - row['Latitude'])**2 + (model["Longitude"] - row["Longitude"])**2)**0.5).min())
# measurements["distance"] = distance
measurements['difference'] = validation

#exclude outlier:






bias = measurements['difference'].mean()
NMB = measurements['difference'].sum() / measurements['Air Pollution Level'].sum()
measurements['error'] = (measurements['difference'] - bias).abs() / measurements['Air Pollution Level'] * 100

# Remove outliers with a difference greater than a set threshold
threshold = 80000 # You can adjust this value based on what you consider an outlier
measurements = measurements[measurements['error'].abs() < threshold]
print(measurements['error'].describe())


per_bias = bias / measurements['Air Pollution Level'].mean() * 100

print(f'percentage bias: {per_bias}')


print(f'NMB: {NMB}')
print(f'Bias corrected normalized error: {measurements['error'].mean()}')




plt.plot(measurements['error'])
plt.xlabel('Index for measurement stations (-)')
plt.ylabel('Bias-Corrected Normalized Error (%)')
plt.title('Bias-Corrected Normalized Error for NO2')
plt.show()


# Scatter plot of PM2.5 levels
plt.figure(figsize=(10, 6))
plt.scatter(measurements["Longitude"], measurements["Latitude"], alpha=0.5, c=measurements["error"], cmap="Reds")
plt.colorbar(label="Percentage (%)")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("Percent difference between model and measurements for NO2")
plt.grid(True)
plt.show()

# Create a map centered around Europe
map_center = [50, 10]  # Approximate center of Europe
europe_map = folium.Map(location=map_center, zoom_start=4)

# Prepare data for heatmap (latitude, longitude, PM2.5 level)
heat_data = list(zip(measurements["Latitude"], measurements["Longitude"], measurements["error"]))

# Add heatmap to the map
HeatMap(heat_data, radius=10, blur=15, max_zoom=10).add_to(europe_map)

# Define the map file path
map_path = "Model_Validation/Map_Representation/NO2_Europe_Map.html"

# Overwrite the existing file with the updated map
europe_map.save(map_path)

print(f"Map updated: {map_path}. Open it in your browser to see the latest heatmap.")


# Bias
print(f"Bias: {bias}")

#error

print(f"Error:{measurements['error'].mean()}")