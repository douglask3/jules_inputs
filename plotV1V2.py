from netCDF4 import Dataset
import numpy as np
import iris
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import iris.plot as iplt
import cartopy.crs as ccrs
from libs import git_info
from pdb import set_trace as browser

# Define paths and parameters
raw_output_dir = "outputs/jules/"
arr_output_dir = "outputs/regridded/"
veg_dir        = ["veg1", "veg2"]
file_names     = "gswp3.fluxes."
years          = range(2000,2001)
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

    return output_file(fname_out, dat_variable, lat_variable, lon_variable)       

def coor_range(mn, mx, ncells):
        diff = float(mx  - mn) / ncells
        return np.arange(mn + diff/2, mx , diff)  

def output_file(fname, dat, lat_variable = None, lon_variable = None):
    if len(dat.shape) == 2: dat = np.reshape(dat, [1, dat.shape[0], dat.shape[1]])
     
    if lat_variable is None: lat_variable = coor_range(-90 , 90 , dat.shape[1])
    if lon_variable is None: lon_variable = coor_range(0, 360, dat.shape[2])
    
    rootgrp = Dataset(fname, "w", format="NETCDF4")    
    
    time_dim = True if dat.shape[0] > 1 else False

    if time_dim: time = rootgrp.createDimension("time", dat.shape[0])
    lat = rootgrp.createDimension("lat", dat.shape[1])
    lon = rootgrp.createDimension("lon", dat.shape[2])

    if time_dim: times = rootgrp.createVariable("time","f8",("time",))

    latitudes  = rootgrp.createVariable("latitude","f4",("lat",))
    longitudes = rootgrp.createVariable("longitude","f4",("lon",))
    
    dims = ("time","lat","lon",) if time_dim else ("lat","lon",)
    var = rootgrp.createVariable(varname, "f4", dims)
    
    latitudes [:] = lat_variable
    longitudes[:] = lon_variable
    if time_dim:
        var[:,:,:] = dat
    else:
        var[:,:] = dat[0,:,:]
    rootgrp.close()

    return dat

def openFile(veg, varname, year, month):
    m = str(month) if month > 9 else '0' + str(month)
    fname_in  = raw_output_dir + veg + '/' + file_names + str(year) + m + '.nc'
    fname_out = arr_output_dir + veg + '_' + file_names + varname + str(year) + m + '.nc'  
    print(fname_in)  
    return open_and_regrid_file(fname_in, fname_out, varname)

def lat_size(lat, dlat, dlon = None): 
    from math import pi, sin, pow
    if dlon is None: dlon = dlat 

    lat1 = (lat + dlat / 2.0) * pi / 180.0
    lat2 = (lat - dlat / 2.0) * pi / 180.0

    size_lat = abs(sin(lat1) - sin(lat2))
    size_lat = (dlat/360) * 2.0 * size_lat * pi * pow(6378137.0, 2.0)
    size_lat = abs(size_lat)
    return size_lat

def global_total(dat):
    lat = coor_range(-90 , 90 , dat.shape[1])
    
    tot = 0
    for i in range(0, dat.shape[1]): 
        ar = lat_size(lat[i], 1.25) 
        tot = tot + sum(sum(dat[:, 100, :])) * ar
    
    return tot

veg = veg_dir[1]
annual_average = openFile(veg, varname, 2000, 1)
annual_average[:,:,:] = 0.0
TS = np.zeros([len(years)*12])
t   = 0
for y in years:
    for m in range(1,12):
        t = t + 1
        dat = openFile(veg_dir[1], varname, y, m)
        annual_average = annual_average + dat
        TS[t] = global_total(dat)

annual_average = annual_average * 60 * 60 * 24 * 30        
fname = arr_output_dir + veg + '_' + varname +  '.nc' 
output_file(fname, annual_average)

plotable = iris.load(fname)

browser()
git = 'repo: ' + git_info.url + '\n' + 'rev:  ' + git_info.rev
crs_latlon = ccrs.PlateCarree()
crs_proj   = ccrs.Robinson()
fig = plt.figure()
ax = plt.axes(projection = crs_proj)
ax.set_extent((-180, 170, -65, 90.0), crs=crs_latlon)
ax.coastlines(linewidth=0.75, color='navy')
ax.gridlines(crs=crs_latlon, linestyle='--')

iplt.contourf(plotable[2], [0.0, 0.5, 1, 1.5, 2, 2.5, 3], colors=('#6666ff', '#acff88', 'b'))
plt.gca().coastlines()

iplt.citation(git)

plt.show()
#plt.savefig(fig_name, bbox_inches='tight')

browser()

