from netCDF4 import Dataset
import numpy as np

import iris
import iris.plot as iplt
import iris.quickplot as qplt

import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import cartopy.crs as ccrs

from os.path import isfile

from libs import git_info
from pdb import set_trace as browser

# Define paths and parameters
raw_output_dir = "outputs/jules/"
arr_output_dir = "outputs/regridded/"
veg_dir        = ["veg1", "veg2"]
file_names     = "gswp3.fluxes."
years          = range(2000,2004)

varnames       = [['gpp_gb']]
limits         = [[[0.0, 0.5, 1, 1.5, 2, 2.5, 3], [-0.01, 0.01]]]

remake_files   = False
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

def coor_range(mn, mx, ncells):
        diff = float(mx  - mn) / ncells
        return np.arange(mn + diff/2, mx , diff)

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

def openFile(veg, varname, year, month):
    m = str(month) if month > 9 else '0' + str(month)
    fname_in  = raw_output_dir + veg + '/' + file_names + str(year) + m + '.nc'
    fname_out = arr_output_dir + veg + '_' + file_names + varname + str(year) + m + '.nc'
    
    print(fname_in)
    if not remake_files and isfile(fname_out):
        nc  = Dataset(fname_out,  "r+", format = "NETCDF4")
        dat = nc.variables[ varname   ][:]
    else:
        dat =   open_and_regrid_file(fname_in, fname_out, varname)
    
    return dat

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
    lat = coor_range(-90.0 , 90.0 , dat.shape[1])

    tot = 0
    for i in range(0, dat.shape[1]):
        toti = np.nansum(dat[:, i, :]) * lat_size(lat[i], 1.25)        
        if not np.isnan(toti): tot = tot + toti

    return tot

def compare_variable(varname, limits):
    varname = varname[0]
    def open_variable(veg): 
        annual_average = openFile(veg, varname, 2000, 1)
        annual_average[:,:,:] = 0.0
        TS = np.zeros([len(years)*12])
        t   = 0
        for y in years:
            for m in range(1,13):
                dat = openFile(veg, varname, y, m)
                annual_average = annual_average + dat
                TS[t] = global_total(dat)
                t = t + 1

        annual_average = annual_average * 60 * 60 * 24 * 365 / t
        fname = arr_output_dir + veg + '_' + varname +  '.nc'
        output_file(fname, varname, annual_average)
        
        TS = TS * 60.0 * 60.0 * 24.0 * (365.0/12.0) /1000000000000.0

        return(fname, TS)
   
    aa1, ts1 = open_variable(veg_dir[0])
    aa2, ts2 = open_variable(veg_dir[1])
    #############
    ## Plot    ##
    #############
    ## setup plot
    git = 'repo: ' + git_info.url + '\n' + 'rev:  ' + git_info.rev
    crs_latlon = ccrs.PlateCarree()
    crs_proj   = ccrs.Robinson()

    ## Plot maps
    def plot_map(fname, cmap, limits, title, ssp):
        if len(fname) == 2: 
            plotable = [iris.load_cube(i) for i in fname]
            plotable = plotable[1] - plotable[0]
        else:
            plotable = iris.load_cube(fname)
        
        ax = fig.add_subplot(2, 2, ssp, projection = crs_proj)
        ax.set_extent((-180, 170, -65, 90.0), crs = crs_latlon)
        ax.coastlines(linewidth = 0.75, color = 'navy')
        ax.gridlines(crs = crs_latlon, linestyle = '--')
               
        qplt.contourf(plotable[0], 10, cmap = cmap)
        plt.gca().coastlines()
        plt.title(title)
    
    fig = plt.figure(figsize=(12, 8))
    cmap =  brewer_cmap = mpl_cm.get_cmap('brewer_YlGn_09')
    dcmap =  brewer_cmap = mpl_cm.get_cmap('brewer_PRGn_11')
    plot_map(aa1, cmap, limits[0], 'Veg1', 1)
    plot_map(aa2, cmap, limits[0], 'Veg2', 2)
    plot_map([aa1, aa2], dcmap, limits[1], 'Difference',     3)

    iplt.citation(git)

    ## line plot
    ax = fig.add_subplot(224)
    t = np.arange(years[0], years[-1] + 1, 1.0/12.0)

    ax.get_xaxis().get_major_formatter().set_scientific(False)
    plt.plot(t, ts1, 'r', t, ts2, 'b--')

    fig_name = 'figs/' + varname + '_veg2-veg1_comparison.pdf'
    plt.savefig(fig_name, bbox_inches='tight')


for v, l in zip(varnames, limits): compare_variable(v, l)
