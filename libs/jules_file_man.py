#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
from libs.grid_funs import coor_range
from os.path import isfile
from libs.min_diff_sequance import *
from pdb import set_trace as browser

__author__ = "Douglas Kelley"
__email__  = "douglas.i.kelley@gmail.com"


class jules_vector2grid(object):
    def __init__(self, fname, varname):
        nc    = Dataset(fname,  "r+", format = "NETCDF4")
        '''dat   = nc.variables[ varname   ][:][0, :, :]'''
        lat   = nc.variables['latitude' ][:]
        lon   = nc.variables['longitude'][:]

        self.lat_variable, self.lat_index = min_diff_sequance(lat, -90, 90 )
        self.lon_variable, self.lon_index = min_diff_sequance(lon, 0  , 360)     
        
    def  open_and_regrid_file(self, fname_in, fname_out, varname, Fun = sum, perLength = True):
        
        nc    = Dataset(fname_in,  "r+", format = "NETCDF4")
        dat   = nc.variables[ varname   ][:]

        n = dat.shape[0]
        if Fun is not None: dat = np.apply_along_axis(Fun, 0, dat)
        if n == 0:
            dn = n
            browser()
        if perLength: dat = dat/n
        if len(dat.shape) == 2: dat = np.reshape(dat, [1, 1, dat.shape[1]])
        
        dat_variable = np.zeros([ dat.shape[0], len(self.lat_variable), len(self.lon_variable)])
        dat_variable[:, :, :]   = np.NAN

        for i in range(1, dat.shape[2]):
            x = self.lon_index[i - 1]; y = self.lat_index[i - 1]
            dat_variable[:, y, x] = dat[:, 0, i]

        return output_file(fname_out, varname, dat_variable, self.lat_variable, self.lon_variable)

def open_file(fname_in, fname_out, varname, remake_files = False, grid = None):    
    print(fname_in)
    if not remake_files and isfile(fname_out):
        nc  = Dataset(fname_out,  "r+", format = "NETCDF4")
        dat = nc.variables[ varname   ][:]
    else:
        if grid is None: grid = jules_vector2grid(fname_in,  varname)
        dat = grid.open_and_regrid_file(fname_in, fname_out, varname)
    return(dat, grid)


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
