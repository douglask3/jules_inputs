import numpy as np

def coor_range(mn, mx, ncells):
        diff = float(mx  - mn) / ncells
        return np.arange(mn + diff/2, mx , diff)


def lat_size(lat, dlat, dlon = None):
    from math import pi, sin, pow
    if dlon is None: dlon = dlat

    lat1 = (lat + dlat / 2.0) * pi / 180.0
    lat2 = (lat - dlat / 2.0) * pi / 180.0

    size_lat = abs(sin(lat1) - sin(lat2))
    size_lat = (dlat/360) * 2.0 * size_lat * pi * pow(6378137.0, 2.0)
    size_lat = abs(size_lat)
    return size_lat


def global_total(dat):
    lat = coor_range(-90.0 , 90.0 , dat.shape[1])

    tot = 0
    for i in range(0, dat.shape[1]):
        toti = np.nansum(dat[:, i, :]) * lat_size(lat[i], 1.25)        
        if not np.isnan(toti): tot = tot + toti

    return tot
