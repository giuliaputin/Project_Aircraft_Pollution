import xarray as xr
import pandas as pd
import folium
import numpy as np
from folium.plugins import HeatMap
import os

# Load the dataset
ds = xr.open_dataset(os.path.join(os.path.dirname(__file__), "..", "..", 'raw_data', 'model', 'PM25.annual.ON.nc4'))

# Extract latitude, longitude, and PM2.5 data
lat = ds['lat'].values
lon = ds['lon'].values
pm25 = ds['PM25'].values

# Convert to Pandas DataFrame
df = pd.DataFrame({
    "Latitude": np.repeat(lat, len(lon)),  # Repeat lat values
    "Longitude": np.tile(lon, len(lat)),   # Tile lon values
    "PM25": pm25.flatten()                 # Flatten PM2.5 values
})

# Remove NaN values
df = df.dropna()

# Display the full data in table form (save to CSV if needed)
df.to_csv("pm25_data.csv", index=False)
print(df.head(20))  # Show the first 20 rows

# Create an interactive map centered over Europe
m = folium.Map(location=[50, 10], zoom_start=4, tiles="CartoDB dark_matter")

# Aggregating data into a heatmap-friendly format
# Average PM2.5 per unique location (lat, lon) to reduce the clutter
df_grouped = df.groupby(['Latitude', 'Longitude'], as_index=False).agg({'PM25': 'mean'})

# Prepare heatmap data
heat_data = df_grouped[['Latitude', 'Longitude', 'PM25']].values.tolist()

# Adjusted HeatMap with a better radius and blur for visibility
HeatMap(heat_data, radius=20, blur=25, min_opacity=0.3).add_to(m)

# Save map as HTML file and open it in a browser
m.save("PM25_Europe_Map.html")
print("Interactive map saved as 'PM25_Europe_Map.html'. Open it in a browser!")
