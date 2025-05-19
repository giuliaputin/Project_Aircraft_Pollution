import xarray as xr
import plotly.express as px
import pandas as pd

# Open the dataset
ds = xr.open_dataset('O3.JUL.OFF.nc4')
print(ds)

# Select the variable for ozone concentration
var = 'SpeciesConc_O3'  # Corrected variable name for ozone
da = ds[var]

# Average over the vertical dimension (lev) to get the total or averaged O3 concentration at each location
da_avg = da.mean(dim='lev')

# Select data at a specific time, using nearest available time point
da_time = da_avg.sel(time='2019-01-15', method='nearest')

# Create a DataFrame from the selected data for Plotly
data = da_time.to_dataframe().reset_index()

# Create an interactive map using Plotly
fig = px.scatter_geo(data, lat='lat', lon='lon', color='SpeciesConc_O3',
                     hover_name='SpeciesConc_O3',
                     color_continuous_scale='Viridis',
                     projection='natural earth',
                     title="Ozone (O₃) Distribution over Europe on 2019-01-15",
                     labels={'SpeciesConc_O3': 'Ozone Concentration'})

# Show the map
fig.show()
