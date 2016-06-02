###############################
## CFG                       ##
###############################
from netCDF4 import Dataset
import numpy as np

file  = 'outputs/gc3_orca1_mask.nc'
input = Dataset(file, "r+", format = "NETCDF4")

###############################
## make lat/lon arrays       ##
###############################
## Open
lon        = input.variables['longitude'][:]
lat = lat0 = input.variables[ 'latitude'][:]

## Convert 2 arrays
def conc(x0, y):   
    x = np.zeros((y.shape[0], x0.shape[0]))
    for i in range(y.shape[0]): x[i, ] = x0
    return(x)       

lat = conc(lat, lon )
lon = conc(lon, lat0)

###############################
## Output                    ##
###############################
def gridVar(x, nm):
    var = input.createVariable(nm, "f4", ("latitude", "longitude"))
    var[:, :] = x

gridVar(lon  , 'grid_lon')
gridVar(lat.T, 'grid_lat')

input.close()

