import pandas as pd
import matplotlib.pyplot as plt
import folium
from folium.plugins import HeatMap
import os

# Define the file path
file_path = "Measurements_Data/PM25.csv"  # Adjust based on your folder structure

# Check if the file exists
if not os.path.exists(file_path):
    raise FileNotFoundError(f"File not found: {file_path}. Check your folder structure!")

# Load the CSV file
df = pd.read_csv(file_path)

# Check if required columns exist
required_columns = {"Longitude", "Latitude", "Air Pollution Level", "Altitude"}
if not required_columns.issubset(df.columns):
    raise ValueError(f"CSV file is missing one or more required columns: {required_columns}")

long_1 = [] # y<0
lat_1 = [] # y<0
pm25_1 = [] # y<0
long_2 = [] # 0<y<50
lat_2 = [] # 0<y<50
pm25_2 = [] # 0<y<50
long_3 = [] # 50<y<100 
lat_3 = [] # 50<y<100 
pm25_3 = [] # 50<y<100 
long_4 = [] # 100<y<250
lat_4 = [] # 100<y<250
pm25_4 = [] # 100<y<250
long_5 = [] # 250<y<1000
lat_5 = [] # 250<y<1000
pm25_5 = [] # 250<y<1000
long_6 = [] # 1000<y
lat_6 = [] # 1000<y
pm25_6 = [] # 1000<y

for i in range(len(df["Altitude"])):
    if df["Altitude"][i]<=0:
        long_1.append(df["Longitude"][i])
        lat_1.append(df["Latitude"][i])
        pm25_1.append(df["Air Pollution Level"][i])
    elif 0<df["Altitude"][i]<=50:
        long_2.append(df["Longitude"][i])
        lat_2.append(df["Latitude"][i])
        pm25_2.append(df["Air Pollution Level"][i])
    elif 50<df["Altitude"][i]<=100:
        long_3.append(df["Longitude"][i])
        lat_3.append(df["Latitude"][i])
        pm25_3.append(df["Air Pollution Level"][i])
    elif 100<df["Altitude"][i]<=250:
        long_4.append(df["Longitude"][i])
        lat_4.append(df["Latitude"][i])
        pm25_4.append(df["Air Pollution Level"][i])
    elif 250<df["Altitude"][i]<=1000:
        long_5.append(df["Longitude"][i])
        lat_5.append(df["Latitude"][i])
        pm25_5.append(df["Air Pollution Level"][i])
    else: 
        long_6.append(df["Longitude"][i])
        lat_6.append(df["Latitude"][i])
        pm25_6.append(df["Air Pollution Level"][i])

range = int(input("Choose altitude range:\n 1: y<0\n 2: 0<y<50\n 3: 50<y<100\n 4: 100<y<250\n 5: 250<y<1000\n 6: 1000<y\n",))
if range==1:
    long_in = long_1
    lat_in = lat_1
    pm_in = pm25_1
if range==2:
    long_in = long_2
    lat_in = lat_2
    pm_in = pm25_2
if range==3:
    long_in = long_3
    lat_in = lat_3
    pm_in = pm25_3
if range==4:
    long_in = long_4
    lat_in = lat_4
    pm_in = pm25_4
if range==5:
    long_in = long_5
    lat_in = lat_5
    pm_in = pm25_5
if range==6:
    long_in = long_6
    lat_in = lat_6
    pm_in = pm25_6
# Scatter plot of PM2.5 levels inputted range 
plt.figure(figsize=(10, 6))
plt.scatter(long_in, pm_in, alpha=0.5, c=pm_in, cmap="Reds")
plt.colorbar(label="PM2.5 Concentration (µg/m³)")
plt.xlabel("Longitude")
plt.ylabel("PM2.5 Concentration (µg/m³)")
plt.title(f"PM2.5 Concentration Across Europe for Altitude Range {range}")
plt.grid(True)
plt.show()

# Create a map centered around Europe
map_center = [50, 10]  # Approximate center of Europe
europe_map = folium.Map(location=map_center, zoom_start=4)

# Prepare data for heatmap (latitude, longitude, PM2.5 level) range 6
heat_data = list(zip(lat_in, long_in, pm_in))

# Add heatmap to the map
HeatMap(heat_data, radius=10, blur=15, max_zoom=10).add_to(europe_map)

# Define the map file path
map_path = "Processing_altitude/Map_Representation/PM25_Altitude_Map.html"

# Overwrite the existing file with the updated map
europe_map.save(map_path)

print(f"Map updated: {map_path}. Open it in your browser to see the latest heatmap.")