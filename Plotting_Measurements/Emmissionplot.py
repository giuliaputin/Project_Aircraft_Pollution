import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import folium
from folium.plugins import HeatMap
import os
import pandas as pd

# Load Masses dataset
masses_ds = xr.open_dataset('raw_data/emissions/AvEmMasses.nc4')
fluxes_ds = xr.open_dataset('raw_data/emissions/AvEmFluxes.nc4')
# Select variable to plot, e.g., NO2
variable = 'NO2'



data = masses_ds[variable].values.flatten()  # Flatten the


# Coordinates explicitly defined
lats = masses_ds['lat']
lons = masses_ds['lon']

lon_grid, lat_grid = np.meshgrid(lons, lats)

lons = lon_grid.flatten()
lats = lat_grid.flatten()


print('This is the data', data)
print('This is the lats', lats)
print('This is the lons', lons)

# Prepare data for heatmap (latitude, longitude, NO₂ level)
heat_data = list(zip(lats, lons, data))

print(heat_data)

# Create a map centered around Europe
map_center = [50, 10]  # Rough center of Europe
europe_map = folium.Map(location=map_center, zoom_start=4)

# Add heatmap to the map
HeatMap(heat_data, radius=10, blur=15, max_zoom=10).add_to(europe_map)



# Add airport locations to the map

airportscsv = "Airports_file/EU_airports.csv"  

# Check if the file exists before reading
if not os.path.exists(airportscsv):
    raise FileNotFoundError(f"File not found: {airportscsv}. Check your folder structure!")

# Load the CSV file
airports = pd.read_csv(airportscsv)

# Check if required columns exist
required_columns = {"latitude_deg", "longitude_deg", "name", "type"}
if not required_columns.issubset(airports.columns):
    raise ValueError(f"CSV file is missing one or more required columns: {required_columns}")

#Grab specific columns from the airports dataframe
airports = airports[["latitude_deg", "longitude_deg", "name", "type"]]


# Add airport markers to the map
for row in airports.itertuples(index=False):
    lat = row.latitude_deg
    lon = row.longitude_deg
    name = row.name
    type = row.type
    
    if type == "large_airport":
        folium.Marker(
            location=[lat, lon],
            popup=name,
            icon=folium.Icon(color="blue", icon="plane", prefix="fa")
        ).add_to(europe_map) 

    #elif type == "medium_airport":
    #    folium.Marker(
    #        location=[lat, lon],
    #        popup=name,
    #        icon=folium.Icon(color="green", icon="plane", prefix="fa")
    #    ).add_to(europe_map)

    #elif type == "small_airport":
    #    folium.Marker(
    #        location=[lat, lon],
    #        popup=name,
    #        icon=folium.Icon(color="red", icon="plane", prefix="fa")
    #    ).add_to(europe_map)  



# Save and display map
map_path = "Plotting_Measurements/Map_Representation/AvEmFluxes.html"
europe_map.save(map_path)
print(f"Map saved as {map_path}. Open it in your browser to view the heatmap.")
