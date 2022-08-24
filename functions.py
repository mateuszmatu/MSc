import os
import re

def create_directory(date):
    if not os.path.exists('flow_maps'):
        os.mkdir('flow_maps')
    path = f'flow_maps/{date}'
    if not os.path.exists(path):
        os.mkdir(path)
    return path

def files_from_thredds(date):
    date_regex = rf'{date}'
    arr = []
    with open('thredds_urls.txt') as file:
        lines = file.readlines()
        for line in lines:
            if re.findall(date_regex, line):
                line = re.sub('\n', '', line)
                arr.append(line)
    arr.sort()
    return arr

def correct_file(arr, member):
    if member >= 0 and member < 6:
        file = arr[0]
        _member = member
    elif member >= 6 and member < 12:
        file = arr[1]
        _member = member - 6
    elif member >= 12 and member < 18:
        file = arr[2]
        _member = member - 12
    elif member >= 18 and member < 24:
        file = arr[3]
        _member = member - 18
    return file, _member

def name_from_lon_lat(lon,lat):
    name = f'lon{str(lon[0])}_{str(lon[1])}_lat{str(lat[0])}_{str(lat[1])}'
    name = re.sub('\.','-',name)
    return name
