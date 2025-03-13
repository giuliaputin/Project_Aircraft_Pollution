import pandas as pd
import matplotlib.pyplot as plt
import folium
from folium.plugins import HeatMap
import os

# Define the file path
file_path = "Measurements_Data/O3.csv"   # Adjust based on your folder structure

# Check if the file exists before reading
if not os.path.exists(file_path):
    raise FileNotFoundError(f"File not found: {file_path}. Check your folder structure!")

# Load the CSV file
df = pd.read_csv(file_path)

# Check if required columns exist
required_columns = {"Longitude", "Latitude", "Altitude", "Air Pollution Level"}
if not required_columns.issubset(df.columns):
    raise ValueError(f"CSV file is missing one or more required columns: {required_columns}")

temp = []
layer = []
final = []

for j in range(10):
    for i in range(len(df["Altitude"])):
        if df["Altitude"][i] > float(50 * (j)) and df["Altitude"][i] < float(50 * (j + 1)):
            temp = [df["Longitude"][i], df["Latitude"][i], df["Air Pollution Level"][i]]
            layer.append(temp)
            temp = []
    final.append(layer)
    layer = []


lon = []
lat = []
quality = []
count = -1

for k in final:
    count += 1
    for l in k:
        lon.append(l[0])
        lat.append(l[1])
        quality.append(l[2])

    map_center = [50, 10]  # Rough center of Europe
    europe_map = folium.Map(location=map_center, zoom_start=4)

    # Prepare data for heatmap (latitude, longitude, O3 level)
    heat_data = list(zip(lat, lon, quality))

    # Add heatmap to the map
    HeatMap(heat_data, radius=10, blur=15, max_zoom=10).add_to(europe_map)

    # Save and display map
    map_path = "Processing_altitudemeasurements/Map_Representation/O3_Altitude_Map" + "_" + str(count) + ".html"
    
    try:
        europe_map.save(map_path)
    except:
        x = open(map_path, "x")
        x.close
        europe_map.save(map_path)

    print(f"Map saved as {map_path}. Open it in your browser to view the heatmap.")

