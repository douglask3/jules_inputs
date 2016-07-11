from netCDF4 import Dataset
import numpy as np
from libs.grid_funs import coor_range
from os.path import isfile
from libs.min_diff_sequance import *

def open_and_regrid_file(fname_in, fname_out, varname, Fun = sum, perLength = True):
    nc    = Dataset(fname_in,  "r+", format = "NETCDF4")
    dat   = nc.variables[ varname   ][:]
    lat   = nc.variables['latitude' ][:]
    lon   = nc.variables['longitude'][:]

    n = dat.shape[0]
    if Fun is not None: dat = np.apply_along_axis(Fun, 0, dat)
    if n == 0:
        dn = n
        browser()
    if perLength: dat = dat/n
    if len(dat.shape) == 2: dat = np.reshape(dat, [1, 1, dat.shape[1]])

    lat_variable, lat_index = min_diff_sequance(lat, -90 , 90 )
    lon_variable, lon_index = min_diff_sequance(lon, 0, 360)
    dat_variable            = np.zeros([ dat.shape[0], len(lat_variable), len(lon_variable)])
    dat_variable[:, :, :]   = np.NAN

    for i in range(1, dat.shape[2]):
        x = lon_index[i - 1]; y = lat_index[i - 1]
        dat_variable[:, y, x] = dat[:, 0, i]

    return output_file(fname_out, varname, dat_variable, lat_variable, lon_variable)

def open_file(fname_in, fname_out, varname, remake_files = False):    
    print(fname_in)
    if not remake_files and isfile(fname_out):
        nc  = Dataset(fname_out,  "r+", format = "NETCDF4")
        dat = nc.variables[ varname   ][:]
    else:
        dat =   open_and_regrid_file(fname_in, fname_out, varname) 
    return dat


def output_file(fname, varname, dat, lat_variable = None, lon_variable = None):
    if len(dat.shape) == 2: dat = np.reshape(dat, [1, dat.shape[0], dat.shape[1]])
    if lat_variable is None: lat_variable = coor_range(-90 , 90 , dat.shape[1])
    if lon_variable is None: lon_variable = coor_range(0, 360, dat.shape[2])

    rootgrp = Dataset(fname, "w", format="NETCDF4")

    time = rootgrp.createDimension("time", dat.shape[0])
    lat  = rootgrp.createDimension("lat", dat.shape[1])
    lon  = rootgrp.createDimension("lon", dat.shape[2])

    times      = rootgrp.createVariable("time","f8",("time",))
    latitudes  = rootgrp.createVariable("lat","f4",("lat",))
    longitudes = rootgrp.createVariable("lon","f4",("lon",))

    longitudes.lon_name         = 'Longitude'
    longitudes.axis             = "X"
    longitudes.standard_name    = "longitude"
    longitudes.units            = "degrees_east"

    latitudes.lon_name          = 'Latitude'
    latitudes.axis              = "Y"
    latitudes.standard_name     = "latitude"
    latitudes.units             = "degrees_north"

    dims = ("time","lat","lon",)
    var = rootgrp.createVariable(varname, "f4", dims)

    latitudes [:] = lat_variable
    longitudes[:] = lon_variable
    var[:,:,:]    = dat
    
    rootgrp.close()

    return dat
