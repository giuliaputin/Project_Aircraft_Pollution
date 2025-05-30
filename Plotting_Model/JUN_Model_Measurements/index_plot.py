# Import libraries
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
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

standard_density = 44.6  # mol/m^3
molar_mass_NO2 = 46.01  # g/mol
molar_mass_O3 = 48  # g/mol

# Calculate total pollution index (as per your formula)
total_index = (
    0.5 * o3_diff.da_diff * standard_density * molar_mass_O3 * 1e6
    + (25 / 60) * no2_diff.da_diff * standard_density * molar_mass_NO2 * 1e6
    + (25 / 15) * pm25_diff.da_diff
)

# Extract lat/lon
lats = datasets['O3']['lat'].values
lons = datasets['O3']['lon'].values
nlat, nlon = total_index.shape

# Set global font size for plots
plt.rcParams.update({'font.size': 40})

# Start plotting with Cartopy
fig = plt.figure(figsize=(12, 8))
ax = plt.axes(projection=ccrs.EqualEarth(central_longitude=10))

# Add base features
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.5, edgecolor='darkgrey')
ax.coastlines(resolution='50m', linewidth=0.5, color='white')

# Plot total pollution index
im = ax.pcolormesh(
    lons, lats, total_index,
    transform=ccrs.PlateCarree(),
    cmap='RdYlGn_r'
)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05)
cbar.ax.tick_params(labelsize=20)

# ========== Add EU + UK + France + Norway boundary line ==========

# Read countries shapefile using Cartopy
shapefile = shpreader.natural_earth(resolution='110m', category='cultural', name='admin_0_countries')
reader = shpreader.Reader(shapefile)
countries = list(reader.records())

# List of ISO3 codes for EU countries + UK + France + Norway
eu_uk_plus = {
    'AUT', 'BEL', 'BGR', 'HRV', 'CYP', 'CZE', 'DNK', 'EST', 'FIN', 'FRA',  # Added FRA here
    'DEU', 'GRC', 'HUN', 'IRL', 'ITA', 'LVA', 'LTU', 'LUX', 'MLT', 'NLD',
    'POL', 'PRT', 'ROU', 'SVK', 'SVN', 'ESP', 'SWE', 'GBR', 'NOR'  # Added NOR here
}

# Extract geometries of selected countries
selected_geoms = [rec.geometry for rec in countries if rec.attributes['ISO_A3'] in eu_uk_plus]

# Plot selected countries boundaries with blue lines (no fill)
ax.add_geometries(
    selected_geoms,
    crs=ccrs.PlateCarree(),
    edgecolor='blue',
    facecolor='none',
    linewidth=2,
    zorder=10
)

# Save the figure
plt.savefig("aircraft_pollution_europe_jul.pdf", format="pdf", bbox_inches='tight')
plt.show()

# --- Folium interactive map creation ---

flat_values = total_index.values.flatten()
lat_grid, lon_grid = np.meshgrid(lats, lons, indexing='ij')
lat_flat = lat_grid.flatten()
lon_flat = lon_grid.flatten()

valid = ~np.isnan(flat_values)
flat_values = flat_values[valid]
lat_flat = lat_flat[valid]
lon_flat = lon_flat[valid]

# Color and opacity functions (same as your code)
def get_color(val, vmin=0.0, vmax=0.6):
    norm = (val - vmin) / (vmax - vmin)
    norm = max(0, min(norm, 1))
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

vmin = np.nanpercentile(total_index.values, 1)
vmax = np.nanpercentile(total_index.values, 99)

lat_res = abs(lats[1] - lats[0])
lon_res = abs(lons[1] - lons[0])

m = folium.Map(location=[50, 10], zoom_start=4, tiles='cartodbpositron')

for lat, lon, val in zip(lat_flat, lon_flat, flat_values):
    color = get_color(val, vmin, vmax)
    opacity = 0.5  # fixed opacity for now

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

# Add airports (make sure path and CSV file are correct)
airportscsv = os.path.join(os.path.dirname(__file__), "..", "..", "Airportdata", "EU_airports.csv")

if not os.path.exists(airportscsv):
    raise FileNotFoundError(f"File not found: {airportscsv}. Check your folder structure!")

airports = pd.read_csv(airportscsv)

required_columns = {"latitude_deg", "longitude_deg", "name", "type"}
if not required_columns.issubset(airports.columns):
    raise ValueError(f"CSV missing columns: {required_columns}")

airports = airports[["latitude_deg", "longitude_deg", "name", "type"]]

# (Optional) Add airport markers - uncomment to enable
# for row in airports.itertuples(index=False):
#     lat = row.latitude_deg
#     lon = row.longitude_deg
#     name = row.name
#     type = row.type
#     icon_color = {'large_airport':'blue', 'medium_airport':'green', 'small_airport':'red'}.get(type, 'gray')
#     folium.Marker(
#         location=[lat, lon],
#         popup=name,
#         icon=folium.Icon(color=icon_color, icon='plane', prefix='fa')
#     ).add_to(m)

# Save Folium map
output_path = "air_quality_map.html"
m.save(output_path)
