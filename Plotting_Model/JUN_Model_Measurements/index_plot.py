# Import libraries
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import folium
import os
import pandas as pd

# Import precomputed difference arrays
import O3_plotting_JUL_diff as o3_diff
import NO2_plotting_JUL_diff as no2_diff
import PM25_plotting_JUL_diff as pm25_diff

# Create dictionary of datasets
datasets = {
    'O3': o3_diff.da_diff,
    'NO2': no2_diff.da_diff,
    'PM25': pm25_diff.da_diff
}

standard_density = 44.6 #mol/m^3
molar_mass_NO2 = 46.01 #g/mol
molar_mass_O3 = 48 #g/mol

#Define scaled index function
def assign_index_scaled(concentration):
    """
    Scales concentration to a 0–100 range based on percentiles,
    improving contrast while avoiding outlier distortion.
    """
    c = concentration.copy()
    min_val = np.nanpercentile(c, 1)
    max_val = np.nanpercentile(c, 99)
    scaled = (c - min_val) / (max_val - min_val)
    scaled = np.clip(scaled, 0, 1)
    return scaled * 100

# # Calculate index for each pollutant
# indexes = {}
# for pol in datasets:
#     indexes[pol] = assign_index_scaled(datasets[pol])
#
# # Combine indexes and scale again to highlight contrasts
# total_index = indexes['O3']*0.5 + indexes['NO2']*(25/60) + indexes['PM25']*(25/15)
# total_index = assign_index_scaled(total_index)

total_index = 0.5 * o3_diff.da_diff * standard_density * molar_mass_O3 * (1000000) + (25/60) *  no2_diff.da_diff * standard_density * molar_mass_NO2 * (1000000) + (25/15) * pm25_diff.da_diff

# Extract lat/lon
lats = datasets['O3']['lat'].values
lons = datasets['O3']['lon'].values
nlat, nlon = total_index.shape

# Find top 10 most polluted grid points
flat = total_index.values.flatten()
top10_idx = np.argpartition(flat, -1)[-1:]
top10_idx = top10_idx[np.argsort(flat[top10_idx])[::-1]]

# Get lat/lon and index values of those points
top10_coords = [(lats[i // nlon], lons[i % nlon], flat[i]) for i in top10_idx]

# Set global font size
plt.rcParams.update({'font.size': 40})  # You can increase this value for even larger text

# Plotting
fig = plt.figure(figsize=(12, 8))
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Plot total pollution index
im = ax.pcolormesh(
    datasets['O3']['lon'], datasets['O3']['lat'], total_index,
    transform=ccrs.PlateCarree(),
    cmap='RdYlGn_r'
)

# Add colorbar and title with larger font sizes
cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05)
cbar.set_label('Total Pollution Index (Aviation Contribution)', fontsize=23)
cbar.ax.tick_params(labelsize=20)

ax.set_title('Overall Aircraft Pollution Impact Over Europe in January', fontsize=25)

# Save figure
# plt.savefig(r"C:\Users\unigi\Downloads\aircraft_pollution_europe.svg", format="svg", bbox_inches='tight')
# plt.savefig("aircraft_pollution_europe_jan.pdf", format="pdf", bbox_inches='tight')

# Blue markers and labels
# for lat, lon, index_value in top10_coords:
#     ax.plot(lon, lat, 'bo', markersize=5, transform=ccrs.PlateCarree())
#     ax.text(lon + 0.5, lat + 0.5, f'{index_value:.1f}', color='blue', fontsize=10,
#             transform=ccrs.PlateCarree())
#     print(lat, lon, index_value)

plt.show()


flat_values = total_index.values.flatten()
lat_grid, lon_grid = np.meshgrid(lats, lons, indexing='ij')
lat_flat = lat_grid.flatten()
lon_flat = lon_grid.flatten()

valid = ~np.isnan(flat_values)
flat_values = flat_values[valid]
lat_flat = lat_flat[valid]
lon_flat = lon_flat[valid]

# Define color and opacity functions
def get_color(val, vmin=0.0, vmax=0.6):
    norm = (val - vmin) / (vmax - vmin)
    norm = max(0, min(norm, 1))

    # Green → Yellow → Red
    if norm < 0.5:
        ratio = norm / 0.5
        r = int(26 + ratio * (255 - 26))
        g = int(150 + ratio * (255 - 150))
        b = int(65 + ratio * (191 - 65))
    else:
        ratio = (norm - 0.5) / 0.5
        r = int(255 + ratio * (215 - 255))
        g = int(255 + ratio * (48 - 255))
        b = int(191 + ratio * (39 - 191))

    return f'#{r:02x}{g:02x}{b:02x}'


def get_opacity(val, vmin, vmax, min_opacity=0.05, max_opacity=0.9):
    norm = (val - vmin) / (vmax - vmin)
    norm = max(0, min(norm, 1))
    return min_opacity + norm * (max_opacity - min_opacity)

# Min/max for scaling
vmin = np.nanpercentile(total_index.values, 1)
vmax = np.nanpercentile(total_index.values, 99)

# Grid resolution
lat_res = abs(lats[1] - lats[0])
lon_res = abs(lons[1] - lons[0])

# Create Folium map
m = folium.Map(location=[50, 10], zoom_start=4, tiles='cartodbpositron')

# Add grid rectangles
for lat, lon, val in zip(lat_flat, lon_flat, flat_values):
    color = get_color(val, vmin, vmax)
    opacity = 0.5   #get_opacity(val, vmin, vmax)

    bounds = [
        [lat - lat_res / 2, lon - lon_res / 2],
        [lat + lat_res / 2, lon + lon_res / 2]
    ]

    folium.Rectangle(
        bounds=bounds,
        color='gray',
        weight=0.5,
        fill=True,
        fill_color=color,
        fill_opacity=opacity,
        tooltip=f"{val:.2f}"
     ).add_to(m)


# Add airport locations to the map
airportscsv = os.path.join(os.path.dirname(__file__), "..", "..", "Airportdata", "EU_airports.csv")

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
    # if type == "large_airport":
    #     folium.Marker(
    #         location=[lat, lon],
    #         popup=name,
    #         icon=folium.Icon(color="blue", icon="plane", prefix="fa")
    #     ).add_to(m) 

    # elif type == "medium_airport":
    #     folium.Marker(
    #         location=[lat, lon],
    #         popup=name,
    #         icon=folium.Icon(color="green", icon="plane", prefix="fa")
    #     ).add_to(m)

    # elif type == "small_airport":
    #     folium.Marker(
    #         location=[lat, lon],
    #         popup=name,
    #         icon=folium.Icon(color="red", icon="plane", prefix="fa")
    #     ).add_to(m)


# Save and display
output_path = "air_quality_map.html"
m.save(output_path)