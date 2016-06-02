cd outputs
./../libs/conv2nc_ltol.tcl ../data/ancil_files/*
./../libs/conv2nc_ltol.tcl ../data/gc3_orca1_mask
cd ..
python addLatLon.py
ncwa -a t outputs/qrparm.veg.nc -O outputs/qrparm.veg.nc
