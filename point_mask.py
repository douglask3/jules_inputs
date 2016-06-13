#################################################
## cfg                                         ##
#################################################
## Libraries
import csv
from netCDF4 import Dataset
import numpy as np
import iris
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import iris.plot as iplt
import cartopy.crs as ccrs
from pdb import set_trace as browser
from libs import git_info

## File info
file_in  = 'docs/points_for_mask2.csv'
file_out = 'outputs/points_file.txt'

file_lsm  = 'outputs/gc3_orca1_mask.nc'
varn_lsm  = 'lsm'

file_lsm_plot = 'data/gc3_orca1_mask'
fig_name  = 'figs/mask_and_points.pdf'

#################################################
## Open                                        ##
#################################################
## Mask
lsm  = Dataset(file_lsm,  "r+", format = "NETCDF4")
lon  = lsm.variables ['longitude'][:]
lat  = lsm.variables [ 'latitude'][:]
lsm  = lsm.variables [  varn_lsm ][:]
mask = iris.load_cube(file_lsm_plot)

## Lat and lons for point mask from file
lonLat = np.genfromtxt(file_in, delimiter = ',')
lonLat[lonLat[:, 0] < 0.0, 0] += 360

#################################################
## Find closest on map                         ##
#################################################
def round2grid(x, y): 
    return([y[np.argmin(abs(y - i))] for i in x])

lonLat[1:, 0] = round2grid(lonLat[1:, 0], lon)
lonLat[1:, 1] = round2grid(lonLat[1:, 1], lat)

#################################################
## Output                                      ##
#################################################
git = 'repo: ' + git_info.url() + '\n' + 'rev:  ' + git_info.rev
## Write to file
np.savetxt(file_out, lonLat[:,[1,0]], fmt = '%.4f', delimiter = '\t')

## Plot
crs_latlon = ccrs.PlateCarree()
crs_proj   = ccrs.Robinson()
fig = plt.figure()
ax = plt.axes(projection = crs_proj)
ax.set_extent((-180, 170, -65, 90.0), crs=crs_latlon)
ax.coastlines(linewidth=0.75, color='navy')
ax.gridlines(crs=crs_latlon, linestyle='--')

iplt.contourf(mask, [-0.5, 0.5, 1.5], colors=('#6666ff', '#acff88', 'b'))
lonLat = crs_proj.transform_points(crs_latlon,lonLat[1:,0],lonLat[1:,1])# lonLat[1:,]

plt.scatter(lonLat[:,0],lonLat[:,1], c = '#aa0000', marker = 'x')
plt.gca().coastlines()

iplt.citation(git)

plt.savefig(fig_name, bbox_inches='tight')

