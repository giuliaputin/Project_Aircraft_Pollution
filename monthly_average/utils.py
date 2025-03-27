import xarray as xr
import os

def differencer(type, month, var):
    dsoff = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." ,  "raw_data", "model", f"{type}.{month}.OFF.nc4") )
    dson = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." , "raw_data", "model", f"{type}.{month}.ON.nc4") )
    daoff = dsoff[var]
    daon = dson[var]

    daSurfoff = daoff.isel(lev=0).mean(dim = "time")
    daSurfon = daon.isel(lev=0).mean(dim = "time")
    daDiff = daSurfon - daSurfoff

    return daDiff

def adder(type, month, var):
    dsoff = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." ,  "raw_data", "model", f"{type}.{month}.OFF.nc4") )
    dson = xr.open_dataset( os.path.join(os.path.dirname(__file__), ".." , "raw_data", "model", f"{type}.{month}.ON.nc4") )
    daoff = dsoff[var]
    daon = dson[var]

    daSurfoff = daoff.isel(lev=0).mean(dim = "time")
    daSurfon = daon.isel(lev=0).mean(dim = "time")
    daDiff = daSurfon - daSurfoff

    return daDiff