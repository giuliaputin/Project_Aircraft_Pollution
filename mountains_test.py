import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import geopandas as gpd
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

# Create a Cartopy map with satellite imagery
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': ccrs.PlateCarree()})

# Use Cartopy's Google Maps (satellite) tile service
tiles = cimgt.GoogleTiles(style='satellite')
ax.add_image(tiles, 8)  # Zoom level 8 for a moderate view

# Add a shapefile (or manually created shapely geometries)
# For this example, we'll create a polygon using Shapely
polygon = Polygon([(0, 30), (0, 50), (60, 50), (60, 90)])
ax.set_extent([0, 40, 29.5, 90], crs=ccrs.PlateCarree())  # Adjust extent for the region

# Use Shapely geometry to plot the region on the map
x, y = polygon.exterior.xy
ax.plot(x, y, color='blue', linewidth=2, marker='o', markersize=4, transform=ccrs.PlateCarree())

# Add a point using Shapely
point = Point(-77, 33)
ax.scatter(point.x, point.y, color='red', transform=ccrs.PlateCarree())

# Add title and show the plot
ax.set_title("Satellite Map with Shapely Overlay")
plt.show()
