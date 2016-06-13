#!/usr/bin/env python

###############################################
## cfg                                       ##
###############################################
## import libs
import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import numpy as np
from netCDF4 import Dataset
from os import listdir, getcwd, mkdir, path
from pdb import set_trace
from libs import git_info

# Define paths and parameters
root_dir = "data/"
clim_dir = "half_deg_files/"
outp_dir = "outputs/"

###############################################
## Load mask and remove coordinate system    ##
###############################################
mask_file = root_dir + "gc3_orca1_mask"
mask_list = iris.load( mask_file)
mask = mask_list[0]

mask.coord("latitude").coord_system=None
mask.coord("longitude").coord_system=None

###############################################
## Re-grid climate files and make drive_file ##
###############################################
drive_date = []
fnames = []

info = ('Regridded to gc3_orca1_mask using ' + git_info.url +  
        ', rev number ' + git_info.rev + 
        ' by Douglas Kelley, douglas.i.kelley@gmail.com')

for input_file in listdir(root_dir + clim_dir):
    if input_file == 'GSWP3.BC.SWdown.3hrMap.2010.N96e.nc': continue
    
    ## Cal filenames
    
    file = root_dir + clim_dir + input_file
    print(file)
   
    ## load cube
    cube = iris.load_cube( file) 
    date = str(cube.coord('time')[0])[10:29]

    output_dir  = outp_dir + cube.var_name + "/"
    output_file = output_dir + "regridded_" + input_file
    
    ## Find new drive_file entry
    if date not in drive_date: 
        drive_date.append(date)
        fname = '%vv/%vv'.join(output_file.split(cube.var_name))
        fnames.append(fname) 

    ## Regrid and output climate
    output_data = cube.regrid( mask, iris.analysis.Nearest())
    
    if not path.isdir(output_dir): mkdir(output_dir)
    iris.save(output_data, output_file)
    
    output = Dataset(output_file, 'r+', format='NETCDF4')
    output.regridding = info
    output.close()

## Output drive_file
drive_date = sorted(drive_date)
fnames     = sorted(fnames)

drive_info = '\n'.join("'" + getcwd() + '/' + a + "', '" + b + "'" for a,b in zip(fnames, drive_date))
drive_file = open("outputs/drive_file.txt", "w")

drive_file.write(drive_info)
drive_file.close()


