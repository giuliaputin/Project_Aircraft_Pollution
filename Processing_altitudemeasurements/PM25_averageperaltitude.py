import matplotlib
matplotlib.use('TkAgg')  # Use 'TkAgg' if 'QtAgg' doesn't work
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import numpy as np
import os

# Define the file path
#file_path = "Measurements_Data/O3.csv"   # Adjust based on your folder structure
file_path = "C:/Users/joaob/OneDrive/Documents/GitHub/Project_Aircraft_Pollution/Measurements_Data/PM25.csv"

# Check if the file exists before reading
if not os.path.exists(file_path):
    raise FileNotFoundError(f"File not found: {file_path}. Check your folder structure!")

# Load the CSV file
df = pd.read_csv(file_path)

# Check if required columns exist
required_columns = {"Longitude", "Latitude", "Air Pollution Level"}
if not required_columns.issubset(df.columns):
    raise ValueError(f"CSV file is missing one or more required columns: {required_columns}")

maxlon = 0.0
maxlat = 0.0

minlon = 0.0
minlat = 0.0

for i in range(len(df["Longitude"])):
    if (df["Longitude"][i]) >= maxlon:
        maxlon = df["Longitude"][i]
    if (df["Latitude"][i]) >= maxlat:
        maxlat = df["Latitude"][i]
    if (df["Longitude"][i]) <= minlon:
        minlon = df["Longitude"][i]
    if (df["Latitude"][i]) <= minlat:
        minlat = df["Latitude"][i]

print(f"Longitude: {minlon} to {maxlon}, Latitude: {minlat} to {maxlat}")


# Load your dataset (ensure it has 'latitude', 'longitude', and 'AQI' columns)
# Example: df = pd.read_csv("air_quality_data.csv")
'''
df = pd.DataFrame({
    'latitude': np.random.uniform(-21.3344, 78.9067, 1000),
    'longitude': np.random.uniform(-61.7272, 55.5826, 1000),
    'AQI': np.random.randint(10, 300, 1000)  # Simulated AQI values
})
'''
df = pd.DataFrame({"Lon":df["Longitude"], "Lat": df["Latitude"], "AQI":df["Air Pollution Level"]})


# Define grid resolution (change as needed)
GRID_SIZE = 1.5 # 1° x 1° grid

# Compute grid indices
df['lat_bin'] = np.round(df['Lat'] / GRID_SIZE) * GRID_SIZE
df['lon_bin'] = np.round(df['Lon'] / GRID_SIZE) * GRID_SIZE

min_lat = 30.0  # Minimum latitude
max_lat = 75.0   # Maximum latitude
min_lon = -15.0  # Minimum longitude
max_lon = 40.0   # Maximum longitude

# Filter the dataframe to keep only the rows within the specified range
df = df[(df['Lat'] >= min_lat) & (df['Lat'] <= max_lat) &
        (df['Lon'] >= min_lon) & (df['Lon'] <= max_lon)]

# Group by grid cells and compute the average AQI
grid_avg = df.groupby(['lat_bin', 'lon_bin'])['AQI'].mean().reset_index()

# Rename columns for clarity
grid_avg.rename(columns={'lat_bin': 'latitude', 'lon_bin': 'longitude', 'AQI': 'avg_AQI'}, inplace=True)

# Display results


# Optional: Save to CSV for visualization
# grid_avg.to_csv("averaged_air_quality.csv", index=False)


# Load dataset (Ensure your dataset has latitude, longitude, and avg_AQI columns)
df = grid_avg.rename(columns={"lat_bin": "latitude", "lon_bin": "longitude", "AQI": "avg_AQI"})
# Load dataset (Ensure your dataset has latitude, longitude, and avg_AQI columns)

# Set up the figure and axis using Cartopy
fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

# Add map features (coastlines, borders, land, ocean)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

# Scatter plot of AQI values
sc = ax.scatter(
    df["longitude"], df["latitude"], c=df["avg_AQI"], cmap="coolwarm",
    edgecolor='k', s=70, alpha=0.75  # Adjust size and transparency
)

# Add a color bar
cbar = plt.colorbar(sc, ax=ax, orientation="vertical", label="Air Quality Index (AQI)")

# Set labels and title
ax.set_title("Average PM2.5 Concentration in The Atmosphere", fontsize=14)
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")

# Show the plot


aqi_grid = grid_avg.pivot_table(index="latitude", columns="longitude", values="avg_AQI")

# Convert to NumPy array
aqi_values = aqi_grid.to_numpy()

# Get latitude & longitude bins
lat_bins = aqi_grid.index.values
lon_bins = aqi_grid.columns.values

# Create a meshgrid
lon_mesh, lat_mesh = np.meshgrid(lon_bins, lat_bins)

# Plot using pcolormesh
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
c = ax.pcolormesh(lon_mesh, lat_mesh, aqi_values, cmap="coolwarm", transform=ccrs.PlateCarree())

# Add colorbar and map features
fig.colorbar(c, ax=ax, label="Average PM2.5 (\u03BCg/m3) ")
ax.coastlines()
ax.gridlines(draw_labels=True)

plt.title("Average PM2.5 Concentration in The Atmosphere")
plt.show()