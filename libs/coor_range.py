def coor_range(mn, mx, ncells):
        diff = float(mx  - mn) / ncells
        return np.arange(mn + diff/2, mx , diff)
