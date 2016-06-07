import csv
from netCDF4 import Dataset
import numpy as np

file_in  = 'docs/points_for_mask.csv'
file_out = 'outputs/points_file.txt'

file_lsm  = 'outputs/gc3_orca1_mask.nc'
varn_lsm  = 'lsm'

lsm  = Dataset(file_lsm,  "r+", format = "NETCDF4")
lon = lsm.variables ['longitude'][:]
lat = lsm.variables [ 'latitude'][:]

lonLat = np.genfromtxt(file_in, delimiter = ',')

lonLat[lonLat[:, 0] < 0.0, 0] += 360

def round2grid(x, y): 
    return([y[np.argmin(abs(y - i))] for i in x])

lonLat[1:, 0] = round2grid(lonLat[1:, 0], lon)
lonLat[1:, 1] = round2grid(lonLat[1:, 1], lat)

np.savetxt(file_out, lonLat[:,[1,0]], fmt = '%.4f', delimiter = '\t')

