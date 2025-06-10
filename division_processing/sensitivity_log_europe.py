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


#FONT
font = {'family' : 'DejaVu Sans',
        'size'   : 22}

matplotlib.rc('font', **font)

def safe_lognorm(data, margin_factor=5):
    vmin = np.nanmin(data)
    vmax = np.nanmax(data)

    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin <= 0:
        raise ValueError("Invalid or non-positive data for LogNorm")

    if vmin == vmax:
        # Expand the range around the single value
        vmin /= margin_factor
        vmax *= margin_factor

    return matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)

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
monthname = ["January", "July"]

types = {
    "Aerosol": {
        "PM25": {
            "emittants": ["NO2"],
            "unit": 1e-3,
            "name": r'PM$_{2.5} \:$'
        },
        "AerMassNIT": {
            "emittants": ["NO2"],
            "unit": 1e-3,  # Example, adjust as needed
            "name": r'Inorganic Nitrate Aerosols'
        },
        "AerMassSO4":{
            "emittants": ["FUELBURN"],
            "unit": 1e-3,
            "name": r"Sulfate Aerosol"
        }
        # 'AerMassNH4': ['NO2'],
        # 'AerMassPOA': ['nvPM', 'HC'],
        # 'AerMassBC': ['nvPM']
    },
    # "Aerosol2": {
    #     "AerMassNIT": ["NO2"]
    # },
    "O3": {
        "SpeciesConc_O3": {
            "emittants": ["NO2", "HC", "CO"],
            "unit": 44.6 * 48 * 1e-3,
            "name": r'O$_3 \:$'
        }
    },
    "NO2": {
        "SpeciesConc_NO2": {
            "emittants": ["NO2"],
            "unit": 44.6 * 46.01 * 1e-3,
            "name": r'NO$_2 \:$'
        }
    },
}


percentile = 80

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

for type_, vars in types.items():

    for var, properties in vars.items():
        emittants = properties["emittants"]
        unit = properties["unit"]
        name = properties["name"]

        for m, month in enumerate(months):
            pollutant = differencer(type_, month, var) * unit * area["AREA"] * 128
            print(np.nanmin(area["AREA"].values))
            sum_emittants = adder(emissions, emittants)
            # print(np.any(emittants==0))

            measure = pollutant / sum_emittants

            fig, ax = plt.subplots(
                1,1,
                figsize=[12, 8],
                subplot_kw={"projection": ccrs.EqualEarth()},
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
                draw_labels=False, linewidth=0.3, color="gray", alpha=0.5, linestyle="--"
            )
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {"size": 16}
            gl.ylabel_style = {"size": 16}

            masked_measure = mask_ocean(measure)

            # Step 2: Take absolute value
            abs_measure = abs(masked_measure)

            # Step 3: Compute cutoff and mask values below percentile
            abs_cutoff = np.nanpercentile(abs_measure.values, percentile)
            topx_abs_measure = abs_measure.where(abs_measure.values >= abs_cutoff, np.nan)

            # Step 4: Use the topx mask on the original (signed) data
            topx_masked_measure = masked_measure.where(~np.isnan(topx_abs_measure), np.nan)

            # Step 5: Split into positives and negatives
            positive_measures_topx = topx_masked_measure.where(topx_masked_measure.values > 0)
            negative_measures_topx = abs(topx_masked_measure.where(topx_masked_measure.values < 0))
            
            # print(min(negative_measures.lat))
            # print(min(masked_measure_europe.lat))

            norm1 = safe_lognorm(positive_measures_topx.values)

            # fig.set_size_inches(12, 6)
            if np.any(~np.isnan(negative_measures_topx.values)):
                norm2 = safe_lognorm(negative_measures_topx.values)
                
                negative_measures_topx.plot(
                    ax=ax,
                    cmap="BuGn",
                    transform=ccrs.PlateCarree(),
                    norm=norm2,
                    cbar_kwargs={
                        "label": "Negative Sensitivity",
                        "aspect": 80,
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
                    "label": "Positive Sensitivity",
                    "aspect": 80,
                    "location": "bottom",
                    "pad": 0.01,
                    # "shrink": 0.5,
                },
            )

            # ax[i].set_title(f"{var}, {month}", fontsize=14)
            # plt.tight_layout()
            ax.set_title(f"{name}, {monthname[m]}", fontsize=22)
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
                ),
                dpi = 300
            )

end = time.time()


print(f"Process run in {end - start} s")
