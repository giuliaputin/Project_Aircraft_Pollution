# Import modules
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import folium
from folium.plugins import HeatMap
import pandas as pd
import numpy as np
import xarray as xr
import os

# Open DataSet and print an overview of it
ds = xr.open_dataset(os.path.join(os.path.dirname(__file__), "..",'raw_data', 'emissions', 'AvEmFluxes.nc4'))


# Select a DataArray
var = 'nvPM'
da = ds[var]

# Select data from ground level
flight_level = 15  # Level index
daSurf = da.isel(lev=flight_level)  # Level selection

# Print dataset info
print(daSurf)
print("Min:", daSurf.min().values, "Max:", daSurf.max().values)
print("Mean:", daSurf.mean().values)

# **Apply scaling factor for better readability**
scale_factor = 1e15  # Convert from kg/m²/s to picograms per m²/s (pg/m²/s)
daSurf_scaled = daSurf * scale_factor

# Compute the new mean after scaling
mean_value = float(daSurf_scaled.mean().values)
print(f'\nMean of {var} is = {mean_value:.3f} pg/m²/s')

# **Adjust plot scale dynamically**
vmin = float(daSurf_scaled.min().values)
vmax = float(daSurf_scaled.max().values)

# Prevent vmax from being too small (adjust as needed)
if vmax < vmin + 1e-3:
    vmax = vmin + 1e-3  # Ensure there's enough contrast

# **Plot the scaled data**
fig = plt.figure(figsize=[12, 6])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Use scaled data to ensure visibility
daSurf_scaled.plot(transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)

# Update title to include flight level
plt.title(f'{var} (scaled to pg/m²/s) - Flight Level {flight_level}')
plt.show()


# Extract lat, lon, and emission values
lats = ds['lat'].values
lons = ds['lon'].values
values = daSurf_scaled.values

# Ensure the data is flattened for folium plotting
lat_grid, lon_grid = np.meshgrid(lats, lons, indexing='ij')
lat_flat = lat_grid.flatten()
lon_flat = lon_grid.flatten()
values_flat = values.flatten()

# Remove NaN values for proper plotting
valid_idx = ~np.isnan(values_flat)
data = pd.DataFrame({
    'lat': lat_flat[valid_idx],
    'lon': lon_flat[valid_idx],
    'value': values_flat[valid_idx]
})

data['normalized_value'] = (data['value'] - vmin) / (vmax - vmin)

# Initialize a folium map
m = folium.Map(location=[0, 0], zoom_start=2, tiles='cartodbpositron')

# Add a heatmap
HeatMap(data[['lat', 'lon', 'normalized_value']].values, min_opacity=0, radius=10, blur=15).add_to(m)


# Add airport locations to the map
airportscsv = os.path.join(os.path.dirname(__file__), "..","Airportdata", "EU_airports.csv")

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
    
    # Uncomment this block to add airport markers to the map
    if type == "large_airport":
        folium.Marker(
            location=[lat, lon],
            popup=name,
            icon=folium.Icon(color="blue", icon="plane", prefix="fa")
        ).add_to(m) 
    '''
    elif type == "medium_airport":
        folium.Marker(
            location=[lat, lon],
            popup=name,
            icon=folium.Icon(color="green", icon="plane", prefix="fa")
        ).add_to(m)

    elif type == "small_airport":
        folium.Marker(
            location=[lat, lon],
            popup=name,
            icon=folium.Icon(color="red", icon="plane", prefix="fa")
        ).add_to(m)  '''
        


# Save map to file
m.save(os.path.join(os.path.dirname(__file__), "..", 'Plotting_Emissions', 'Map_Representation', 'plotting_fluxes_PM_NO2_HC.html'))

# Display map (if running in Jupyter Notebook, use `m` to show inline)
print("Map has been saved as 'plotting_fluxes_PM_NO2_HC.html'. Open this file in a browser to view it.")