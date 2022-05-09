import os
from ftle import FTLE
import re
import numpy as np

def generate_LCS(path):
    list_of_files = []

    for root, dirs, files in os.walk(path):
        for file in files:
            list_of_files.append(os.path.join(root,file))
    
    list_of_files.sort()
    regex1 = r'(?<=flow_maps/)(.*?)(?=/)'
    path1 = str(re.findall(regex1, path)[0])
    regex2 = rf'(?<=flow_maps/{path1}/)(.*?)(?=$)'
    path2 = str(re.findall(regex2, path)[0])

    vals = f'flow_maps/{path1}/initial_values.npy'

    if not os.path.exists('LCS_files'):
        os.mkdir('LCS_files')
    spath = f'LCS_files/{path1}'
    if not os.path.exists(spath):
        os.mkdir(spath)
    spath = f'LCS_files/{path1}/{path2}'
    if not os.path.exists(spath):
        os.mkdir(spath)

    for i, file in enumerate(list_of_files):
        ftle = FTLE(file, vals)
        np.save(f'{spath}/ftle_m{i+1}.npy', ftle)

if __name__ == '__main__':
    generate_LCS('flow_maps/20220420/b')