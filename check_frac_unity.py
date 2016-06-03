import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import numpy as np

from pdb import set_trace

file = 'outputs/qrparm.veg.frac.nc'
acceptable_limit = 0.05

dat = iris.load_cube(file)
tot = dat.collapsed('pseudo', iris.analysis.SUM)


limits = [0, 1 - acceptable_limit, 1 + acceptable_limit, 2]
if ((tot.data > limits[3]).any) or ((tot.data < limits[1]).any) :
    ## plot discrepancies
    iplt.contourf(tot, limits)
    plt.gca().coastlines()
    plt.colorbar()
    # Add the contour line levels to the colorbar

report = []
for i in range(tot.data.shape[0]):
    for j in range(tot.data.shape[1]):
        if tot.data[i, j] != 1.0:
            report.append([i, j, tot.data[i,j]])
            dat.data[:, i, j] = dat.data[:, i, j] / tot.data[i, j]

iris.save(dat, file)         
set_trace()

#dat.data = dat.data / tot.data




