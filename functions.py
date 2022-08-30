import os
import re

def create_directory(date):
    '''
        Creates a directory under flow_maps/date/, where date is an argument
    Args:
        date    [str]   :   the name of the directory, usually a date
    Returns:
        path    [str]   :   path to directory
    '''
    if not os.path.exists('flow_maps'):
        os.mkdir('flow_maps')
    path = f'flow_maps/{date}'
    if not os.path.exists(path):
        os.mkdir(path)
    return path

def files_from_thredds(date):
    '''
        Find the threeds files from thredds_urls.txt from the specified date
    Args:
        date    [str]   :   the date from which files should be found
    Returns:
        arr     [list]  :   a list of the file names from the specified date
    '''
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
    '''
        Finds the correct file in the list of files from files_from_thredds(), based on which member is chosen
    Args:
        arr     [list]  :   a list of four filenames
        member  [int]   :   the member to look at
    Returns:
        file    [str]   :   a string of the correct filename
        _member [int]   :   an integer corresponding to the requested member, but in the selected file
    '''
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
    '''
        Creates a string name created from given longitude and latitude
    Args:
        lon     [list]  :   a list of min and max lon
        lat     [list]  :   a list of min and max lat
    Returns:
        name    [str]   :   a string with lons and lats
    '''
    name = f'lon{str(lon[0])}_{str(lon[1])}_lat{str(lat[0])}_{str(lat[1])}'
    name = re.sub('\.','-',name)
    return name
