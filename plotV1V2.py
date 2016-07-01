from netCDF4 import Dataset
import numpy as np
import iris
import matplotlib.pyplot as plt
from pdb import set_trace as browser

# Define paths and parameters
output_dir = "outputs/jules/"
veg_dir    = ["veg1", "veg2"]
file_names = "gswp3.fluxes."
years      = range(2000,2010)
varname    = 'gpp_gb'

###############################################
## Open stuff                                ##
###############################################

files = [output_dir + veg_dir[0] + '/' + file_names + str(i) + '.nc' for i in years]

def min_diff_sequance(x, min_v, max_v):
    d =  np.diff(np.unique(x)).min()
    y = np.arange(min_v + d / 2, max_v + d / 2, d)
    i = (x - min_v - d/2)/d
    i = i.astype(int)[0,:]
    return(y, i)

def openFile(veg, varname, year, month):
    m = str(month) if month > 9 else '0' + str(month)
    
    fname = output_dir + veg + '/' + file_names + str(year) + m + '.nc'
    nc    = Dataset(fname,  "r+", format = "NETCDF4")
    dat   = nc.variables[ varname   ][:]
    lat   = nc.variables['latitude' ][:]
    lon   = nc.variables['longitude'][:]
    
    
    lat_variable, lat_index = min_diff_sequance(lat, -90 , 90 )
    lon_variable, lon_index = min_diff_sequance(lon, 0, 360)
    dat_variable = np.zeros([len(lon_variable), len(lat_variable), dat.shape[0]])
    
    for i in range(1, dat.shape[2]):
        print(i)
        #browser()
        x = lon_index[i - 1]; y = lat_index[i - 1] 
        dat_variable[x, y, :] = dat[:, 0, i]

    # Now either need to write to nc file or add smoothing operation
    # smoothing might be fastest way to get plots for this project

openFile(veg_dir[1], varname, 2000, 1)



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
