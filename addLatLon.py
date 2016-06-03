###############################
## CFG                       ##
###############################
from netCDF4 import Dataset
import numpy as np
from pdb import set_trace

file_lsm  = 'outputs/gc3_orca1_mask.nc'
file_frac = 'outputs/qrparm.veg.frac.nc'

varn_lsm  = 'lsm'
varn_frac = 'field1391'


###############################
## make lat/lon arrays       ##
###############################
## Open
input_lsm  = Dataset(file_lsm,  "r+", format = "NETCDF4")
input_frac = Dataset(file_frac, "r" , format = "NETCDF4")

lon        = input_lsm.variables ['longitude'][:]
lat = lat0 = input_lsm.variables [ 'latitude'][:]
lsm        = input_lsm.variables [ varn_lsm  ][:]
frac       = input_frac.variables[ varn_frac ][:]

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
    var = input_lsm.createVariable(nm, "f4", ("latitude", "longitude"))
    var[:, :] = x

gridVar(lon  , 'grid_lon')
gridVar(lat.T, 'grid_lat')

#lsm  = lsm[0, 0, ]
lsm[0, 0, frac.mask[1,]] = 0.0
gridVar(lsm  , 'grid_lsm')



input_lsm.close()

