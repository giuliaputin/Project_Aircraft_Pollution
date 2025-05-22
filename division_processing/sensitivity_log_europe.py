import matplotlib.pyplot as plt
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import os
import numpy as np
from utils import adder, differencer
import time
import shapely.geometry as sgeom
import cartopy.io.shapereader as shpreader
from shapely.prepared import prep
from matplotlib.patches import Rectangle

# Load the 'admin_0_countries' shapefile (countries of the world)
shpfilename = shpreader.natural_earth(
    resolution="10m", category="cultural", name="admin_0_countries"
)
reader = shpreader.Reader(shpfilename)
records = reader.records()

europe_countries = [
    "Austria",
    "Belgium",
    "Bosnia and Herzegovina",
    "Bulgaria",
    "Croatia",
    "Czech Republic",
    "Denmark",
    "Estonia",
    "Finland",
    "France",
    "Germany",
    "Greece",
    "Hungary",
    "Iceland",
    "Ireland",
    "Italy",
    "Latvia",
    "Lithuania",
    "Luxembourg",
    "Malta",
    "Moldova",
    "Netherlands",
    "North Macedonia",
    "Norway",
    "Poland",
    "Portugal",
    "Romania",
    "Serbia",
    "Slovakia",
    "Slovenia",
    "Spain",
    "Sweden",
    "Switzerland",
    "Ukraine",
    "United Kingdom",
    "Albania",
    "Montenegro",
    "Kosovo",
    "Belarus",
]

# Collect geometries for only the European countries
polygons = []
for record in records:
    if record.attributes["NAME_LONG"] in europe_countries:
        geom = record.geometry
        if isinstance(geom, sgeom.MultiPolygon):
            polygons.extend(geom.geoms)
        else:
            polygons.append(geom)

# Create a MultiPolygon for Europe
europe_geom = sgeom.MultiPolygon(polygons)
prep_europe = prep(europe_geom)

# Rough bounding box for Europe
EUROPE_BOUNDS = {"lon_min": -25, "lon_max": 40, "lat_min": 35, "lat_max": 68}


# ----------------------------------------------------
# Function to mask ocean values from a DataArray
# ----------------------------------------------------
def mask_ocean(dataarray):
    lons, lats = np.meshgrid(dataarray["lon"], dataarray["lat"])
    points = np.vstack([lons.ravel(), lats.ravel()]).T
    mask = np.array([prep_europe.contains(sgeom.Point(x, y)) for x, y in points])
    mask = mask.reshape(lats.shape)
    return dataarray.where(mask)


# ----------------------------------------------------
# Plot only over Europe
# ----------------------------------------------------
start = time.time()
months = ["JAN", "JUL"]

types = {
    "Aerosol": {
        "PM25": ["nvPM"],
        # 'AerMassNIT': ['NO2'],
        # 'AerMassNH4': ['NO2'],
        # 'AerMassPOA': ['nvPM', 'HC'],
        # 'AerMassBC': ['nvPM']
    },
    "O3": {"SpeciesConc_O3": ["NO2", "HC", "CO"]},
    "NO2": {"SpeciesConc_NO2": ["NO2"]},
}


percentile = 0

units = [1e-3, 44.6 * 48 * 1e-3, 44.6 * 46.01 * 1e-3]

emissions = xr.open_dataset(
    os.path.join(
        os.path.dirname(__file__), "..", "raw_data", "emissions", "AvEmMasses.nc4"
    )
)
area = xr.open_dataset(
    os.path.join(
        os.path.dirname(__file__), "..", "raw_data", "model", "daAREA_05x0625_EU.nc4"
    )
)

for i, values in enumerate(types.items()):
    type_, vars = values

    for var, emittants in vars.items():

        for month in months:
            pollutant = differencer(type_, month, var) * units[i] * area["AREA"] * 128
            sum_emittants = adder(emissions, emittants)

            measure = pollutant / sum_emittants

            fig, ax = plt.subplots(
                1,1,
                figsize=[12, 8],
                subplot_kw={"projection": ccrs.PlateCarree()},
                constrained_layout=True,
            )

            # for i, month in enumerate(months):
            ax.add_feature(
                cfeature.BORDERS.with_scale("10m"), linewidth=0.5, edgecolor="darkgrey"
            )
            ax.coastlines(resolution="10m", linewidth=0.5, color="black")

            # Set extent to Europe
            ax.set_extent(
                [
                    EUROPE_BOUNDS["lon_min"],
                    EUROPE_BOUNDS["lon_max"],
                    EUROPE_BOUNDS["lat_min"],
                    EUROPE_BOUNDS["lat_max"],
                ],
                crs=ccrs.PlateCarree(),
            )

            # Add bounding box on map
            europe_rect = Rectangle(
                (EUROPE_BOUNDS["lon_min"], EUROPE_BOUNDS["lat_min"]),
                EUROPE_BOUNDS["lon_max"] - EUROPE_BOUNDS["lon_min"],
                EUROPE_BOUNDS["lat_max"] - EUROPE_BOUNDS["lat_min"],
                linewidth=1.5,
                edgecolor="blue",
                facecolor="none",
                linestyle="--",
                transform=ccrs.PlateCarree(),
            )
            # ax[i].add_patch(europe_rect)

            # Optional: Add gridlines with labels
            gl = ax.gridlines(
                draw_labels=True, linewidth=0.3, color="gray", alpha=0.5, linestyle="--"
            )
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {"size": 9}
            gl.ylabel_style = {"size": 9}

            masked_measure = mask_ocean(measure)
            negative_measures = abs(masked_measure.where(masked_measure.values < 0))
            positive_measures = masked_measure.where(masked_measure.values > 0)

            positive_measures_cutoff = np.nanpercentile(positive_measures.values, percentile)
            positive_measures_topx = positive_measures.where(
                positive_measures.values >= positive_measures_cutoff, np.nan
            )

            negative_measures_cutoff = np.nanpercentile(negative_measures.values, percentile)
            negative_measures_topx = negative_measures.where(
                negative_measures.values >= negative_measures_cutoff, np.nan
            )

            # print(min(negative_measures.lat))
            # print(min(masked_measure_europe.lat))

            vmin = np.nanmin(positive_measures_topx.values)
            vmax = np.nanmax(positive_measures_topx.values)
            norm1 = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)

            # fig.set_size_inches(12, 6)
            if np.any(~np.isnan(negative_measures_topx.values)):
                norm2 = matplotlib.colors.LogNorm(
                    vmin=np.nanmin(negative_measures_topx.values),
                    vmax=np.nanmax(negative_measures_topx.values),
                )
                print(np.nanmin(negative_measures_topx.values), np.nanmax(negative_measures_topx.values))
                negative_measures_topx.plot(
                    ax=ax,
                    cmap="BuGn",
                    transform=ccrs.PlateCarree(),
                    norm=norm2,
                    cbar_kwargs={
                        "label": "Negative values (log scale)",
                        # "aspect": 40,
                        "location": "bottom",
                        "pad": 0.01,
                        # "shrink": 0.5,
                    },
                )

            positive_measures_topx.plot(
                ax=ax,
                transform=ccrs.PlateCarree(),
                norm=norm1,
                cmap="YlOrRd",
                cbar_kwargs={
                    "label": "Positive value (log scale)",
                    # "aspect": 40,
                    "location": "bottom",
                    "pad": 0.01,
                    # "shrink": 0.5,
                },
            )

            # ax[i].set_title(f"{var}, monthly average {month}", fontsize=14)
            # plt.tight_layout()
            ax.set_title("")
            if percentile > 0:
                ifpercen = "_top" + str(100 - percentile)
            else:
                ifpercen = ""
            
            plt.savefig(
                os.path.join(
                    os.path.dirname(__file__),
                    "..",
                    "division_processing",
                    "division_figures",
                    "analysis",
                    f"{type_}_{var}_{month}_timeaveraged_logscale{ifpercen}.png",
                )
            )

end = time.time()
plt.show()

print(f"Process run in {end - start} s")
