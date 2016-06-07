import csv
import pandas as pd

file_in  = 'docs/points_for_mask.csv'
file_out = 'outputs/points_file.txt'

lonLat = pd.read_csv(file_in, sep = ',')
lonLat[[1,0]].to_csv(file_out, header = False, index = False)
