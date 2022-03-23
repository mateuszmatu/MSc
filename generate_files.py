from time import time
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
from netCDF4 import Dataset
import numpy as np
import re
from datetime import timedelta
import os

def initiate_files(grid_size, lon, lat, time_step, file, outfile_name = None):
    outfile_end = ['f', 'b']
    h=2
    if outfile_name == None:
        regex = r'(?<=_eps_eps/barents_eps_)(.*?)(?=T)'
        outfile_name = str(re.findall(regex, file)[0])

    if not os.path.exists('flow_maps'):
        os.mkdir('flow_maps')
    path = f'flow_maps/{outfile_name}'
    if not os.path.exists(path):
        os.mkdir(path)
    if not os.path.exists(f'{path}/f'):
        os.mkdir(f'{path}/f')
    if not os.path.exists(f'{path}/b'):
        os.mkdir(f'{path}/b')

    lons = np.linspace(lon[0], lon[1], grid_size)
    lats = np.linspace(lat[0], lat[1], grid_size)
    lons, lats = np.meshgrid(lons, lats)
    lons = lons.ravel()
    lats = lats.ravel()

    for i in range(2):
        for j in range(6):
            o = OceanDrift(loglevel=20)
            o.set_config('drift:vertical_mixing', False)
            r = reader_netCDF_CF_generic.Reader(file, ensemble_member=j)
            o.add_reader(r)
            o.seed_elements(lons, lats, time=r.start_time+timedelta(hours=h),
                            radius=0, number=grid_size**2)

            if i == 0:
                o.run(duration=timedelta(hours=h), time_step=time_step, time_step_output=3600*h, outfile=f'{path}/f/{outfile_name}_{outfile_end[i]}_m{j}.nc')

            if i == 1:
                o.run(duration=timedelta(hours=h), time_step=-time_step, time_step_output=-3600*h, outfile=f'{path}/b/{outfile_name}_{outfile_end[i]}_m{j}.nc')

if __name__ == '__main__':
    lon = [8,15]
    lat = [70,71]
    url = 'https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_eps/barents_eps_20220323T00Z.nc'
    time_step = 15*60
    grid_size = 10
    initiate_files(grid_size, lon, lat, time_step, url)
    

