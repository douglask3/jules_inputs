#!/usr/bin/env python

## import libs
import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import numpy as np
from os import listdir
import pdb

# Define paths and parameters
root_dir = "data/"
clim_dir = "half_deg_files/"
outp_dir = "outputs/"

show_plots = False

## Load mask and remove coordinate system
mask_file = root_dir + "gc3_orca1_mask"
mask_list = iris.load( mask_file)
mask = mask_list[0]

print( "Original coord_system")
print( mask.coord("latitude").coord_system)
print( mask.coord("longitude").coord_system)

#mask.coord("longitude").units="m"
#mask.coord("latitude").units="m"
mask.coord("latitude").coord_system=None
mask.coord("longitude").coord_system=None

print( "New coord_system")
print( mask.coord("latitude").coord_system)
print( mask.coord("longitude").coord_system)

drive_date = []
fnames = []
for input_file in listdir(root_dir + clim_dir):

    output_file = outp_dir + "regridded_" + input_file
    file = root_dir + clim_dir + input_file
    print(file)
   
    cube = iris.load_cube( file)
    
    date = str(cube.coord('time')[1])[10:29]
    if date not in drive_date: 
        drive_date.append(date)
        fname = '%vv'.join(output_file.split(cube.var_name))
        fnames.append(fname)

    input_data = cube[0]
    output_data = input_data.regrid( mask, iris.analysis.Nearest())
    iris.save(output_data, output_file)


drive_info = '\n'.join("'" + a + "', '" + b + "'" for a,b in zip(fnames, drive_date))

drive_file = open("outputs/drive_file.txt", "w")
drive_file.write(drive_info)
drive_file.close()

#pdb.set_trace()
