from   netCDF4 import Dataset
import numpy   as np

import iris
import iris.plot      as iplt
import iris.quickplot as qplt

import matplotlib.pyplot as plt
import matplotlib.colors as colours
from   pylab             import *
import matplotlib.cm     as mpl_cm
import cartopy.crs       as ccrs

from libs                   import git_info
from libs.grid_funs         import * #coor_range, lat_size, global_total
from libs.jules_file_man    import * #open_and_regrid_file, output_file


from matplotlib.mlab import bivariate_normal

from pdb import set_trace as browser

# Define paths and parameters
raw_output_dir = "outputs/jules3/"
arr_output_dir = "outputs/regridded/"
veg_dir        = ["veg1", "veg2"]
file_names     = "gswp3.fluxes."
years          = range(2000,2010)

varnames       = [['gpp_gb'], ['cs_gb'], ['lit_c_mean'], ['npp_gb'], ['resp_l'], ['resp_w'],
                  ['resp_p_gb'], ['resp_r'], ['resp_s'], ['resp_s_gb']]

limits         = [[[0.0, 0.5, 1, 1.5, 2, 2.5, 3, 4], [-1.5, -1, -0.5, -0.1, 0.1, 0.5, 1, 1.5]], #gpp_b
                  [[0.0, 1, 5, 10, 20, 30, 40], [-20, -16, -8, -4, -1, 1, 4, 8, 16, 20]],        #cs_gb
                  [[0.0, 0.1, 0.5, 1, 1.5, 2, 2.5, 3], [-3, -2, -1, 1, 2, 3]],                  #lit_c_mean
                  [[-20, -15, -10, -5, -2, 1, 0.1, 0.2, 0.5,  1, 1.5, 2], [-0.01, 0.01]],       #npp_gb
                  [[0, 0.1, 0.5, 1, 2, 5, 10, 20, 50],
                        [-50, -20, -10, -5, -1, -0.1, 0.1, 1, 5, 10, 20, 50]],                  #resp_l
                  [[0.0, 0.01, 0.05, 1, 0.1, 0.2, 0.3, 0.4],
                        [-50, -20, -10, -5, -1, -0.1, 0.1, 1, 5, 10, 20, 50]],                  #resp_w
                  [[0.0, 0.01, 0.05, 1, 0.1, 0.2, 0.3, 0.4],
                        [-0.16, -0.12, -0.08, -0.04, 0.04, 0.08, 0.12, 0.16]],                  #resp_p_gb
                  [[0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4], 
                        [-0.16, -0.12, -0.08, -0.04, 0.04, 0.08, 0.12, 0.16]],                  #resp_r 
                  [[0.0, 0.5, 1, 1.5, 2, 2.5, 3], [-0.01, 0.01]],                               #resp_s
                  [[0.0, 0.5, 1, 1.5, 2, 2.5, 3], [-0.01, 0.01]]]                               #resp_s_gb

sec2year       = 60 * 60 * 24 * 365
scaling        = [ sec2year, 1    , 1    , sec2year, sec2year, sec2year, sec2year, sec2year, sec2year]
units          = ['PgC/yr', 'PgC', 'PgC', 'PgC/yr', 'PgC/yr', 'PgC/yr',  'PgC/yr', 'PgC/yr', 'PgC/yr']

remake_files   = False
diagnose_lims  = True
###############################################
## Open stuff                                ##
###############################################
def loadFile(veg, varname, year, month, **kwargs):
    m = str(month) if month > 9 else '0' + str(month)
    fname_in  = raw_output_dir + veg + '/' + file_names + str(year) + m + '.nc'
    fname_out = arr_output_dir + veg + '_' + file_names + varname + str(year) + m + '.nc'
    return open_file(fname_in, fname_out, varname, remake_files, **kwargs)


def to_precision(x,p):
    """
    returns a string representation of x formatted with a precision of p

    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """

    x = float(x)

    if x == 0.:
        return 0.0

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)

    if n < math.pow(10, p - 1):
        e = e -1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= math.pow(10,p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return float("".join(out))

def hist_limits(dat, nlims, symmetrical = True):
    lims = np.percentile(dat.data[~np.isnan(dat.data)], range(0, 100, 100/nlims))

    lims = [to_precision(i, 2) for i in lims]
    lims = unique(lims)
    if symmetrical and min(lims) < 0.0 and max(lims) > 0.0: 
        if abs(min(lims)) > max(lims):
            lims = [i for i in lims if i < 0.0]
        else:
            lims = [i for i in lims if i > 0.0]
        
        lims = sort([lims + [-i for i in lims]])[0]
    if len(lims) == 1: 
        lims = lims[0]
        if lims == 0.0:
            lims = [0.0, 0.001, 0.01]
        else:
            lims = [lims - 0.1 * abs(lims), lims, lims + 0.1 *abs(lims)]
            if lims[0] < 0.0 & lims[1] > 0.0: lims[0] = 0.0
            if lims[2] > 0.0 & lims[1] < 0.0: lims[2] = 0.0
    if len(lims) < 3: lims = [2 * lims[0] - lims[1], lims[0], lims[1], 2* lims[1] - lims[0]]
    return lims

def reverse_colourmap(cmap, name = 'my_cmap_r'):
    """
    In: 
    cmap, name 
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
    """        
    reverse = []
    k = []   

    for key in cmap._segmentdata:    
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:                    
            data.append((1-t[0],t[2],t[1]))            
        reverse.append(sorted(data))    

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL) 
    return my_cmap_r


class MidpointNormalize(colours.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colours.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def compare_variable(varname, limits, scaling, units):
    varname = varname[0]
    def open_variable(veg): 
        annual_average, grid = loadFile(veg, varname, 2000, 1)
            
        annual_average[:,:,:] = 0.0
        TS = np.zeros([len(years)*12])
        t   = 0
        for y in years:
            for m in range(1,13):
                dat, nn = loadFile(veg, varname, y, m, grid = grid)
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
    
    nz = iris.load_cube(aa1).shape[0]
    if  nz== 5:
        figsize = (28, 12)
        fontsize = 16
        px = 3
        py = 5
        labs = [['BL'],['NL'],['C3'],['C4'],['Shrub']]
        ndiagnose_lims = 6
    elif nz == 4:
        figsize = (20, 12)
        fontsize = 16
        px = 3
        py = 4
        labs =[['DPM'], ['RPM'], ['BIO'], ['HUM']]
        ndiagnose_lims = 6
    elif nz == 1:
        figsize = (12, 8)
        fontsize = 20
        px = 2
        py = 2    
        labs = [['']]
        ndiagnose_lims = 6
    else:
        browser()
    cmap1 = mpl_cm.get_cmap('brewer_YlGn_09')
    cmap2 = mpl_cm.get_cmap('brewer_RdYlBu_11') #revserse
    #cmap2 = reverse_colourmap(cmap2)
    dcmap = mpl_cm.get_cmap('brewer_PRGn_11')
    ## Plot maps
    def plot_map(fname , limits, title, ssp, extend = 'max'):
        if len(fname) == 2: 
            plotable = [iris.load_cube(i) for i in fname]
            plotable = plotable[1] - plotable[0]
        else:
            plotable = iris.load_cube(fname)
       
        nz = iris.load_cube(aa1).shape[0]
        ssp = (ssp - 1) * nz 
        for i in range(0, nz):
            ssp =  ssp + 1
            print(ssp)

            ax = fig.add_subplot(px, py, ssp, projection = crs_proj)
            ax.set_extent((-180, 170, -65, 90.0), crs = crs_latlon)
            ax.coastlines(linewidth = 0.5, color = 'navy')
            ax.gridlines(crs = crs_latlon, linestyle = '--')
            
            if diagnose_lims: limits = hist_limits(plotable[i], ndiagnose_lims)   
            if limits[0] < 0.0:
                if limits[-1] <= 0.0:
                    cmap = cmap2 
                    #norm = colours.BoundaryNorm(boundaries = np.append(limits, 0.01), ncolors = 8)
                    norm = colours.BoundaryNorm(boundaries = limits, ncolors = 7)
                else:
                    cmap = dcmap
                    norm = colours.BoundaryNorm(boundaries = limits, ncolors = 12)
                    #norm = MidpointNormalize(midpoint=0.)
            else:
                cmap = cmap1   
                norm = colours.BoundaryNorm(boundaries = limits, ncolors = 8)
            plt.gca().coastlines()
            cbi = qplt.contourf(plotable[i], limits, cmap = cmap, extend = extend, norm = norm) 
            #plt.colorbar(cbi)
            titlei = title + labs[i][0]  
            plt.title(titlei)
        
    fig = plt.figure(figsize = figsize)
    fig.suptitle(varname, fontsize = fontsize)
    
    
    plot_map(aa1, limits[0], 'Veg1', 1)
    plot_map(aa2, limits[0], 'Veg2', 2)
    plot_map([aa1, aa2],  limits[1], 'Difference',     3, 'both')

    iplt.citation(git)

    ## line plot
    if nz == 1:
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
    plt.savefig(fig_name, bbox_inches = 'tight')


for v, l, s, u in zip(varnames, limits, scaling, units): compare_variable(v, l, s, u)
