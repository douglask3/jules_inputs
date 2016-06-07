import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import numpy as np
from pdb import set_trace

# Define paths and parameters
file = "jules/outputs/gswp3.fluxes.2000.nc"

###############################################
## Load mask and remove coordinate system    ##
###############################################
dats = iris.load( file)

set_trace()
