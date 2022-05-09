from click import open_file
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
import numpy as np
import re
from datetime import timedelta
import os
def initiate_files_one_member(lon, lat, time_step, dxdy, duration, date, member, outfile_name = None):
    outfile_end = ['f', 'b']
    member -= 1

    if outfile_name is None:
        outfile_name = date

    if not os.path.exists('flow_maps'):
        os.mkdir('flow_maps')
    path = f'flow_maps/{outfile_name}'
    if not os.path.exists(path):
        os.mkdir(path)
    if not os.path.exists(f'{path}/f'):
        os.mkdir(f'{path}/f')
    if not os.path.exists(f'{path}/b'):
        os.mkdir(f'{path}/b')
    
    date_regex = rf'{date}'
    thredds_files = []
    with open('thredds_urls.txt') as file:
        lines = file.readlines()
        for line in lines:
            if re.findall(date_regex, line):
                line = re.sub('\n', '', line)
                thredds_files.append(line)
    thredds_files.sort()

    if member >= 0 and member < 6:
        file = thredds_files[0]
        _member = member
    elif member >= 6 and member < 12:
        file = thredds_files[1]
        _member = member - 6
    elif member >= 12 and member < 18:
        file = thredds_files[2]
        _member = member - 12
    elif member >= 18 and member < 24:
        file = thredds_files[3]
        _member = member - 18
    
    o = OceanDrift(loglevel=20)
    r = reader_netCDF_CF_generic.Reader(file, ensemble_member = _member)
    o.add_reader(r)
    x = np.arange(lon[0], lon[1], dxdy)
    y = np.arange(lat[0], lat[1], dxdy)
    #proj = r.proj
    lons, lats = np.meshgrid(x,y)
    #lons, lats = proj(X,Y, inverse=True)
    o.seed_elements(lons.ravel(), lats.ravel(), time=r.start_time)
    o.run(duration=timedelta(seconds=duration), time_step=time_step, outfile=f'{path}/b/particle_drift.nc')

if __name__ == '__main__':
    lon = [12.5,16]
    lat = [68,70]

    time_step = 7200
    duration = 60*60*24*20
    dxdy = 0.5
    #initiate_file_all_members(lon, lat, time_step, dxdy, duration, '20220420')
    initiate_files_one_member(lon, lat, time_step, dxdy, duration, '20220419', 1)
    