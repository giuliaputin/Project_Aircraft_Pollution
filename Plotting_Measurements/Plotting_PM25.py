import pandas as pd
import matplotlib.pyplot as plt
import folium
from folium.plugins import HeatMap
import os

# Define the file path
file_path = "raw_data/measurements/PM25.csv"  # Adjust based on your folder structure

# Check if the file exists
if not os.path.exists(file_path):
    raise FileNotFoundError(f"File not found: {file_path}. Check your folder structure!")

# Load the CSV file
df = pd.read_csv(file_path)

# Check if required columns exist
required_columns = {"Longitude", "Latitude", "Air Pollution Level"}
if not required_columns.issubset(df.columns):
    raise ValueError(f"CSV file is missing one or more required columns: {required_columns}")

# Scatter plot of PM2.5 levels
plt.figure(figsize=(10, 6))
plt.scatter(df["Longitude"], df["Air Pollution Level"], alpha=0.5, c=df["Air Pollution Level"], cmap="Reds")
plt.colorbar(label="PM2.5 Concentration (µg/m³)")
plt.xlabel("Longitude")
plt.ylabel("PM2.5 Concentration (µg/m³)")
plt.title("PM2.5 Concentration Across Europe")
plt.grid(True)
plt.show()

# Create a map centered around Europe
map_center = [50, 10]  # Approximate center of Europe
europe_map = folium.Map(location=map_center, zoom_start=4)

# Prepare data for heatmap (latitude, longitude, PM2.5 level)
heat_data = list(zip(df["Latitude"], df["Longitude"], df["Air Pollution Level"]))

# Add heatmap to the map
HeatMap(heat_data, radius=10, blur=15, max_zoom=10).add_to(europe_map)

# Define the map file path
map_path = "Plotting_Measurements/Map_Representation/PM25_Europe_Map.html"

# Overwrite the existing file with the updated map
europe_map.save(map_path)

print(f"Map updated: {map_path}. Open it in your browser to see the latest heatmap.")

