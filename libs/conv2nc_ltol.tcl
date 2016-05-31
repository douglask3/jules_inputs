#!/users/global/rjel/bin/convsh

#  Convsh script conv2nc_1to1.tcl
#
#  Convert each file, into a corresponding netCDF file.
#  File names are taken from the command arguments.
#  For example to convert all PP files in the current directory
#  into netCDF files use the following command:
#
#      ./conv2nc_1to1.tcl *.pp

#  Write out netCDF file
set outformat netcdf

#  Automatically work out input file type
set filetype 0

#  Convert all fields in input files to netCDF
set fieldlist -1

#  Read in each of the input files and write output file

foreach infile $argv {

#  Replace input file extension with .nc to get output filename
   set outfile [file tail [file rootname $infile].nc]

#  Read input file
   readfile $filetype $infile

#  Write out all input fields to a netCDF file
   writefile $outformat $outfile $fieldlist

#  Remove input file information from Convsh's memory
   clearall
}
