from csv import reader
from click import open_file
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
import numpy as np
import re
from datetime import timedelta
import os
import xarray as xr

class initiate_particle_advection:
    def __init__(self, lon, lat, time_step, dxdy, duration, date, hours):
        self.lon = lon
        self.lat = lat
        self.time_step = time_step
        self.dxdy = dxdy
        self.duration = duration
        self.date = date

        x = np.arange(lon[0], lon[1], dxdy)
        y = np.arange(lat[0], lat[1], dxdy)
        self.lons, self.lats = np.meshgrid(x,y)

        if not os.path.exists('flow_maps'):
            os.mkdir('flow_maps')
        self.path = f'flow_maps/{date}'
        if not os.path.exists(self.path):
            os.mkdir(self.path)

        date_regex = rf'{date}'
        self.thredds_files = []
        with open('thredds_urls.txt') as file:
            lines = file.readlines()
            for line in lines:
                if re.findall(date_regex, line):
                    line = re.sub('\n', '', line)
                    self.thredds_files.append(line)
        self.thredds_files.sort()

        self.time_per_sim = self.duration/self.time_step
        self.number_of_sims = int(self.time_per_sim/hours)
    
    def set_up_dictionary(self, number_of_sims):
        #Creating a dictionary with values to be saved
        values = {'duration': self.duration, 'time_step': self.time_step, 'dxdy': self.dxdy, 'initial_lon': self.lons, 'initial_lat': self.lats, 'time': {}}
        for i in range(number_of_sims):
            values['time'][i+1] = {}
            values['time'][i+1]['f'] = {}
            values['time'][i+1]['b'] = {}
        return values

    @property
    def initiate_advection_all_members(self):
        #Run for all members for the current day
        for i in range(1,25):
            self.inititate_advection_one_member(i)

    def inititate_advection_one_member(self, member):
        #Run for a single member of the current day
        member -= 1

        if member >= 0 and member < 6:
            file = self.thredds_files[0]
            _member = member
        elif member >= 6 and member < 12:
            file = self.thredds_files[1]
            _member = member - 6
        elif member >= 12 and member < 18:
            file = self.thredds_files[2]
            _member = member - 12
        elif member >= 18 and member < 24:
            file = self.thredds_files[3]
            _member = member - 18
        
        values = self.set_up_dictionary(self.number_of_sims)
        for i in range(self.number_of_sims):
            for j in range(2):
                o = OceanDrift(loglevel=20)
                r = reader_netCDF_CF_generic.Reader(file, ensemble_member=_member)
                o.add_reader(r)
                if j == 0:
                    o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(seconds=i*self.duration))
                    o.run(duration=timedelta(seconds=self.duration), time_step=self.time_step, time_step_output=self.duration, outfile=f'{self.path}/m{member+1}_f.nc')
                if j == 1:
                    o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(seconds=i*self.duration)+timedelta(seconds=self.duration))
                    o.run(duration=timedelta(seconds=self.duration), time_step=-self.time_step, time_step_output=-self.duration, outfile=f'{self.path}/m{member+1}_b.nc')
            _f = xr.open_dataset(f'{self.path}/m{member+1}_f.nc')
            _b = xr.open_dataset(f'{self.path}/m{member+1}_b.nc')
            os.remove(f'{self.path}/m{member+1}_f.nc')
            os.remove(f'{self.path}/m{member+1}_b.nc')
            values['time'][i+1]['f']['lon'] = np.array(_f['lon'][:,-1])
            values['time'][i+1]['f']['lat'] = np.array(_f['lat'][:,-1])
            values['time'][i+1]['b']['lon'] = np.array(_b['lon'][:,-1])
            values['time'][i+1]['b']['lat'] = np.array(_b['lat'][:,-1])
        np.save(f'{self.path}/values_m{member+1}.npy', values)
        
if __name__ == '__main__':
    lon = [12.5,16]
    lat = [68,70]

    time_step = 900
    duration = 7200
    dxdy = 0.02
    #initiate_file_all_members(lon, lat, time_step, dxdy, duration, '20220420')
    r = initiate_particle_advection(lon,lat,time_step,dxdy,duration,'20220319', hours=4)

    #r.initiate_advection_all_members
    r.inititate_advection_one_member(1)
    r.inititate_advection_one_member(2)

#https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_eps_files/barents_eps_20220420T18Z.nc
#https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_eps_files/barents_eps_20220420T12Z.nc
#https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_eps_files/barents_eps_20220420T06Z.nc
#https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_eps_files/barents_eps_20220420T00Z.nc
    """
    @property
    def _initiate_advection_all_members(self):
        ensemble = 1
        values = self.set_up_dictionary(self.number_of_sims)
        
        for file in self.thredds_files:
            for m in range(6):
                
                for t in range(self.number_of_sims):
                    for i in range(2):
                        o = OceanDrift(loglevel=20)
                        r = reader_netCDF_CF_generic.Reader(file, ensemble_member=m)
                        o.add_reader(r)
                        if i == 0:
                            o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(seconds=t*self.duration))
                            o.run(duration=timedelta(seconds=self.duration), time_step=self.time_step, time_step_output=self.duration, outfile=f'{self.path}/f.nc')
                        if i == 1:
                            o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(seconds=t*self.duration)+timedelta(seconds=self.duration))
                            o.run(duration=timedelta(seconds=self.duration), time_step=-self.time_step, time_step_output=-self.duration, outfile=f'{self.path}/b.nc')
                    _f = xr.open_dataset(f'{self.path}/f.nc')
                    _b = xr.open_dataset(f'{self.path}/b.nc')
                    os.remove(f'{self.path}/f.nc')
                    os.remove(f'{self.path}/b.nc')
                    values['time'][t+1]['f']['lon'] = np.array(_f['lon'][:,-1])
                    values['time'][t+1]['f']['lat'] = np.array(_f['lat'][:,-1])
                    values['time'][t+1]['b']['lon'] = np.array(_b['lon'][:,-1])
                    values['time'][t+1]['b']['lat'] = np.array(_b['lat'][:,-1])
                np.save(f'{self.path}/values_m{ensemble}.npy', values)
                values = self.set_up_dictionary(self.number_of_sims)
                ensemble+=1
    """