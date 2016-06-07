#!/usr/bin/env python

###############################################
## cfg                                       ##
###############################################
## import libs
import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import numpy as np
from os import listdir, getcwd
from pdb import set_trace

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
for input_file in listdir(root_dir + clim_dir):
    if input_file == 'GSWP3.BC.SWdown.3hrMap.2010.N96e.nc': continue
    
    ## Cal filenames
    output_file = outp_dir + "regridded_" + input_file
    file = root_dir + clim_dir + input_file
    print(file)
   
    ## load cube
    cube = iris.load_cube( file) 
    date = str(cube.coord('time')[0])[10:29]
    
    ## Find new drive_file entry
    if date not in drive_date: 
        drive_date.append(date)
        fname = '%vv'.join(output_file.split(cube.var_name))
        fnames.append(fname) 

    ## Regrid and ouput climate
    output_data = cube.regrid( mask, iris.analysis.Nearest())
    iris.save(output_data, output_file)

## Output drive_file
drive_date = sorted(drive_date)
fnames     = sorted(fnames)

drive_info = '\n'.join("'" + getcwd() + '/' + a + "', '" + b + "'" for a,b in zip(fnames, drive_date))
drive_file = open("outputs/drive_file.txt", "w")

drive_file.write(drive_info)
drive_file.close()


