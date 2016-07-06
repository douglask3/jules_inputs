from netCDF4 import Dataset
import numpy as np
import iris
import matplotlib.pyplot as plt
from pdb import set_trace as browser

# Define paths and parameters
raw_output_dir = "outputs/jules/"
arr_output_dir = "outputs/regridded/"
veg_dir        = ["veg1", "veg2"]
file_names     = "gswp3.fluxes."
years          = range(2000,2010)
varname        = 'gpp_gb'

###############################################
## Open stuff                                ##
###############################################


def min_diff_sequance(x, min_v, max_v):
    d =  np.diff(np.unique(x)).min()
    y = np.arange(min_v + d / 2, max_v + d / 2, d)
    i = (x - min_v - d/2)/d
    i = i.astype(int)[0,:]
    return(y, i)

def open_and_regrid_file(fname_in, fname_out, varname, Fun = sum, perLength = True):
    nc    = Dataset(fname_in,  "r+", format = "NETCDF4")
    dat   = nc.variables[ varname   ][:]
    lat   = nc.variables['latitude' ][:]
    lon   = nc.variables['longitude'][:]
   
    n = dat.shape[0]
    if Fun is not None: dat = np.apply_along_axis(Fun, 0, dat)
    if perLength: dat = dat/n
    if len(dat.shape) == 2: dat = np.reshape(dat, [1, 1, dat.shape[1]])

    lat_variable, lat_index = min_diff_sequance(lat, -90 , 90 )
    lon_variable, lon_index = min_diff_sequance(lon, 0, 360)
    dat_variable = np.zeros([ dat.shape[0], len(lat_variable), len(lon_variable)])
        
    for i in range(1, dat.shape[2]):
        x = lon_index[i - 1]; y = lat_index[i - 1] 
        dat_variable[:, y, x] = dat[:, 0, i]
       
    rootgrp = Dataset(fname_out, "w", format="NETCDF4")    
    time_dim = True if dat_variable.shape[0] > 1 else False

    if time_dim: time = rootgrp.createDimension("time", dat_variable.shape[0])
    lat = rootgrp.createDimension("lat", dat_variable.shape[1])
    lon = rootgrp.createDimension("lon", dat_variable.shape[2])

    if time_dim: times = rootgrp.createVariable("time","f8",("time",))

    latitudes  = rootgrp.createVariable("latitude","f4",("lat",))
    longitudes = rootgrp.createVariable("longitude","f4",("lon",))
    
    dims = ("time","lat","lon",) if time_dim else ("lat","lon",)
    var = rootgrp.createVariable(varname, "f4", dims)

    latitudes [:] = lat_variable
    longitudes[:] = lon_variable
    if time_dim:
        var[:,:,:] = dat_variable
    else:
        var[:,:] = dat_variable[0,:,:]
    rootgrp.close()
    return(dat_variable)


def openFile(veg, varname, year, month):
    m = str(month) if month > 9 else '0' + str(month)
    fname_in  = raw_output_dir + veg + '/' + file_names + str(year) + m + '.nc'
    fname_out = arr_output_dir + veg + '_' + file_names + str(year) + m + '.nc'  
    print(fname_in)  
    return open_and_regrid_file(fname_in, fname_out, varname)


dat = openFile(veg_dir[1], varname, 2000, 1)
nt  = len(years)*12
dat = np.zeros([nt, dat.shape[1], dat.shape[2]])
t   = 0
for y in years:
    for m in range(1,12):
        t = t + 1
        dat[t, :, : ] = openFile(veg_dir[1], varname, y, m)

browser()

#def dailySum(x):
#    browser()

#def openFile(veg, varname, year):
#    fname = output_dir + veg + '/' + file_names + str(year) + '.nc'
#    dat =  Dataset(fname,  "r+", format = "NETCDF4")
#    dat =  dat.variables[varname][:]
#    dat = np.apply_along_axis(dailySum, 0, dat).shape 
#    browser()

#def openVeg(veg, varname):
#    for y in years: openFile(veg, varname, y)


#veg1 = openVeg(veg_dir[0], varname)
