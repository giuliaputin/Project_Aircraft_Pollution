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
ds = xr.open_dataset(os.path.join('raw_data', 'emissions', 'AvEmFluxes.nc4'))
print(ds)

# Select a DataArray
var = 'HC'
da = ds[var]

# Select data from ground level at a specific time
daSurf = da.isel(lev=0) #.sel(time='2019-01-15')

# **Apply scaling factor for better readability**
scale_factor = 1e15  # Convert from kg/m²/s to picograms per m²/s (pg/m²/s)
daSurf_scaled = daSurf * scale_factor

# Print the mean value of this selection
print(f'\nMean of {var} is = {float(daSurf.mean().values):.1f} {da.attrs["units"]}')

print("Dataset Structure:\n", ds)
print("\nSelected Variable:", var)
print("\nAttributes:\n", da.attrs)

print("\nDataArray Summary:\n", daSurf)
print("\nData Values:\n", daSurf.values)

print("\nMin:", daSurf.min().values, "Max:", daSurf.max().values)
print("\nMean (Direct):", daSurf.mean().values)
print("\nSum of All Values:", daSurf.sum().values)

# Check NaN count
print("\nNumber of NaNs:", daSurf.isnull().sum().values)


# Plot figure using matplotlib and cartopy. The use of some cosmetic options is illustrated here
fig = plt.figure(figsize=[12, 6])
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')
daSurf.plot(transform=ccrs.PlateCarree(), vmin=0, vmax=50)
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

# Initialize a folium map
m = folium.Map(location=[0, 0], zoom_start=2, tiles='cartodbpositron')

# Add a heatmap
HeatMap(data[['lat', 'lon', 'value']].values, min_opacity=0.2, radius=10, blur=15).add_to(m)


# Add airport locations to the map
airportscsv = "Airportdata/EU_airports.csv"

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
    
    '''	# Uncomment this block to add airport markers to the map
    if type == "large_airport":
        folium.Marker(
            location=[lat, lon],
            popup=name,
            icon=folium.Icon(color="blue", icon="plane", prefix="fa")
        ).add_to(m) 

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
        ).add_to(m)  '
        '''


# Save map to file
m.save('Plotting_Emissions/Map_Representation/plotting_fluxes_template.html')

# Display map (if running in Jupyter Notebook, use `m` to show inline)
print("Map has been saved as 'plotting_fluxes_template.html'. Open this file in a browser to view it.")