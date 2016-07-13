from   netCDF4 import Dataset
import numpy   as np

import iris
import iris.plot      as iplt
import iris.quickplot as qplt

import matplotlib.pyplot as plt
from   pylab             import *
import matplotlib.cm     as mpl_cm
import cartopy.crs       as ccrs

from libs                   import git_info
from libs.grid_funs         import * #coor_range, lat_size, global_total
from libs.jules_file_man    import * #open_and_regrid_file, output_file

from pdb import set_trace as browser

# Define paths and parameters
raw_output_dir = "outputs/jules3/"
arr_output_dir = "outputs/regridded/"
veg_dir        = ["veg1", "veg2"]
file_names     = "gswp3.fluxes."
years          = range(2000,2010)

varnames       = [['resp_r'], ['gpp_gb'], ['cs_gb'], ['lit_c_mean'], ['npp_gb'], ['resp_l'], 
                  ['resp_p_gb'], ['resp_s'], ['resp_s_gb']]

limits         = [[[0.0, 0.5, 1, 1.5, 2, 2.5, 3, 4], [-1.5, -1, -0.5, -0.1, 0.1, 0.5, 1, 1.5]], #gpp_b
                  [[0.0, 1, 5, 10, 20, 30, 40], [-20, -16, -8, -4, -1, 1, 4, 8, 16, 20]],        #cs_gb
                  [[0.0, 0.1, 0.5, 1, 1.5, 2, 2.5, 3], [-3, -2, -1, 1, 2, 3]],                  #lit_c_mean
                  [[-20, -15, -10, -5, -2, 1, 0.1, 0.2, 0.5,  1, 1.5, 2], [-0.01, 0.01]],       #npp_gb
                  [[0, 0.1, 0.5, 1, 2, 5, 10, 20, 50],
                        [-50, -20, -10, -5, -1, -0.1, 0.1, 1, 5, 10, 20, 50]],                  #resp_l
                  [[0.0, 0.01, 0.05, 1, 0.1, 0.2, 0.3, 0.4],
                        [-0.16, -0.12, -0.08, -0.04, 0.04, 0.08, 0.12, 0.16]],                  #resp_p_gb
                  [[0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4], 
                        [-0.16, -0.12, -0.08, -0.04, 0.04, 0.08, 0.12, 0.16]],                  #resp_r 
                  [[0.0, 0.5, 1, 1.5, 2, 2.5, 3], [-0.01, 0.01]],                               #resp_s
                  [[0.0, 0.5, 1, 1.5, 2, 2.5, 3], [-0.01, 0.01]]]                               #resp_s_gb

sec2year       = 60 * 60 * 24 * 365
scaling        = [sec2year, 1    , 1    , sec2year, sec2year, sec2year, sec2year, sec2year]
units          = ['PgC/yr', 'PgC', 'PgC', 'PgC/yr', 'PgC/yr', 'PgC/yr', 'PgC/yr', 'PgC/yr']

remake_files   = False
diagnose_lims  = True
###############################################
## Open stuff                                ##
###############################################
def loadFile(veg, varname, year, month):
    m = str(month) if month > 9 else '0' + str(month)
    fname_in  = raw_output_dir + veg + '/' + file_names + str(year) + m + '.nc'
    fname_out = arr_output_dir + veg + '_' + file_names + varname + str(year) + m + '.nc'
    return open_file(fname_in, fname_out, varname, remake_files)

def compare_variable(varname, limits, scaling, units):
    varname = varname[0]
    def open_variable(veg): 
        annual_average = loadFile(veg, varname, 2000, 1)
            
        annual_average[:,:,:] = 0.0
        TS = np.zeros([len(years)*12])
        t   = 0
        for y in years:
            for m in range(1,13):
                dat = loadFile(veg, varname, y, m)
                annual_average = annual_average + dat
                TS[t] = global_total(dat)
                t = t + 1

        annual_average = annual_average * scaling / t
        fname = arr_output_dir + veg + '_' + varname +  '.nc'
        output_file(fname, varname, annual_average)
        
        TS = TS * scaling /(12 * 1000000000000.0)

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
    def plot_map(fname, cmap, limits, title, ssp, extend = 'max'):
        if len(fname) == 2: 
            plotable = [iris.load_cube(i) for i in fname]
            plotable = plotable[1] - plotable[0]
        else:
            plotable = iris.load_cube(fname)
       
        browser() 
        ax = fig.add_subplot(2, 2, ssp, projection = crs_proj)
        ax.set_extent((-180, 170, -65, 90.0), crs = crs_latlon)
        ax.coastlines(linewidth = 0.75, color = 'navy')
        ax.gridlines(crs = crs_latlon, linestyle = '--')
        
        if diagnose_lims: limits = 10     
        plt.gca().coastlines() 
        qplt.contourf(plotable[0], limits, cmap = cmap, extend = extend)       
        plt.title(title)
    
    fig = plt.figure(figsize=(12, 8))
    fig.suptitle(varname, fontsize=20)
    
    cmap =  brewer_cmap = mpl_cm.get_cmap('brewer_YlGn_09')
    dcmap =  brewer_cmap = mpl_cm.get_cmap('brewer_PRGn_11')
    plot_map(aa1, cmap, limits[0], 'Veg1', 1)
    plot_map(aa2, cmap, limits[0], 'Veg2', 2)
    plot_map([aa1, aa2], dcmap, limits[1], 'Difference',     3, 'both')

    iplt.citation(git)

    ## line plot
    ax = fig.add_subplot(224)
    t = np.arange(years[0], years[-1] + 1, 1.0/12.0)
    
    ax.get_xaxis().get_major_formatter().set_scientific(False)

    plt.plot(t, ts1, 'r', label="Veg 1")
    plt.plot(t, ts2, 'b--', label="Veg 2")
    #plt.plot(t, ts1, 'r', t, ts2, 'b--')
    plt.ylabel(units)

    plt.legend(bbox_to_anchor=(0., 0.898, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)

    fig_name = 'figs/' + varname + '_veg2-veg1_comparison.pdf'
    plt.savefig(fig_name, bbox_inches='tight')


for v, l, s, u in zip(varnames, limits, scaling, units): compare_variable(v, l, s, u)
