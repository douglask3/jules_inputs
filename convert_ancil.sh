cd outputs
./../libs/conv2nc_ltol.tcl ../data/ancil_files/*
./../libs/conv2nc_ltol.tcl ../data/gc3_orca1_mask
cd ..
python addLatLon.py

ncwa -a t outputs/qrparm.veg.frac.nc -O outputs/qrparm.veg.frac.nc
ncwa -a t outputs/qrparm.soil.nc -O outputs/qrparm.soil.nc
ncwa -a surface outputs/qrparm.soil.nc -O outputs/qrparm.soil.nc

python check_frac_unity.py
