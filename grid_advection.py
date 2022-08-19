from csv import reader
from tracemalloc import start
from click import open_file
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
import numpy as np
import re
from datetime import timedelta
import os
import xarray as xr
import functions as f
from datetime import datetime


class particle_grid_displacement:
    def __init__(self, lon, lat, time_step, sep, dur, date, hours):
        '''
            Class for placing particles in a grid and then letting them drift with the flow backwards and/or forwards in time.
        Args:
            lon         [list]    : list of min and max longitude
            lat         [list]    : list of min and max latitude
            time_step   [int]     : time step in seconds
            sep         [float]   : initial separation between particles
            dur         [int]     : duration for how long particles should drift
            date        [str]     : the date for the dataset
            hours       [int]     : hours after dataset start 
        '''

        self.lon = lon
        self.lat = lat
        self.ts = time_step
        self.sep = sep
        self.dur = dur
        self.date = date
        self.h = hours

        x = np.arange(lon[0], lon[1], self.sep)
        y = np.arange(lat[0], lat[1], self.sep)
        self.lons, self.lats = np.meshgrid(x, y)

        self.path = f.create_directory(date)
        self.tf = f.files_from_thredds(date)
        self.name = f.name_from_lon_lat(lon, lat)

    @property
    def _set_up_dictionary(self):
        '''
            Creates a dictionary in which values are saved
        Returns:
            values [dict] : a dictionary containing values to be saved and empty spots to be filled out later
        '''

        values = {'duration': self.dur, 'time_step': self.ts, 'in_separation': self.sep, 'in_lons': self.lons, 'in_lats': self.lats, 'f': {}, 'b': {}}
        return values
    
    def displacement_one_member(self, member, outfile = None, forwards = True, backwards = True):
        '''
            Seeds particles in a grid based on specifications, and lets them drift backwards and 
            forwards with the flow for a specified duration.
        Args:
            member      [int]    : the eps member to be used for particle advection
            outfile     [str]    : a string with the outfile name, default None
            forwards    [bool]   : True if particles should be advected forwards in time, False otherwise
            backwards   [bool]   : True if particles should be advected backwards in time, False otherwise
        Returns:
            string      [str]    : a string containing the path to file and filename
        '''

        corr = f.correct_file(self.tf, member)
        file = corr[0]
        _member = corr[1]

        values = self._set_up_dictionary
        if outfile is None:
            string = f'{self.path}/{self.name}_h{self.h}_m{member}.npy'
        elif outfile is not None:
            string = f'{outfile}.npy'

        for i in range(2):
            o = OceanDrift(loglevel=20)
            r = reader_netCDF_CF_generic.Reader(file, ensemble_member=_member)
            o.add_reader(r)
            if i == 0 and forwards is True:
                o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(hours=self.h))
                o.run(duration=timedelta(hours=self.dur), time_step=timedelta(seconds=self.ts), time_step_output=timedelta(hours=self.dur), outfile=f'{self.path}/m{member}_f.nc')
            if i == 1 and backwards is True:
                o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(hours=self.h))
                o.run(duration=timedelta(hours=self.dur), time_step=timedelta(seconds=-self.ts), time_step_output=timedelta(hours=-self.dur), outfile=f'{self.path}/m{member}_b.nc')
       
        _f = xr.open_dataset(f'{self.path}/m{member}_f.nc')
        _b = xr.open_dataset(f'{self.path}/m{member}_b.nc')
        os.remove(f'{self.path}/m{member}_f.nc')
        os.remove(f'{self.path}/m{member}_b.nc')

        if forwards is True:
            values['f']['lon'] = np.array(_f['lon'][:,-1])
            values['f']['lat'] = np.array(_f['lat'][:,-1])
            if backwards is False:
                values['b']['lon'] = np.zeros_like(_f['lon'][:,-1])
                values['b']['lat'] = np.zeros_like(_f['lat'][:,-1])

        if backwards is True:
            values['b']['lon'] = np.array(_b['lon'][:,-1])
            values['b']['lat'] = np.array(_b['lat'][:,-1])
            if forwards is False:
                values['f']['lon'] = np.zeros_like(_b['lon'][:,-1])
                values['f']['lat'] = np.zeros_like(_b['lat'][:,-1])
        
        np.save(string, values)

        return string

    @property
    def displacement_all_members(self, num_of_members = 24):
        '''
            Runs the displacement_one_member() function for all members in the eps.
        Args:
            num_of_members   [int] : number of members in the eps. Default 24 ---> number of members in the barents 2.5 eps.
        '''
        for i in range(num_of_members):
            self.displacement_one_member(i)
        
#****************************************************************************
#
# >>>>> TO BE CLEANED LATER <<<<<<
#
#****************************************************************************
class initiate_particle_advection:
    def __init__(self, lon = None, lat = None, time_step = None, dxdy = 0.1, duration = None, date = None, hours = 2, start_time = 0):
        self.lon = lon
        self.lat = lat
        self.time_step = time_step
        self.dxdy = dxdy
        self.duration = duration
        self.date = date
        self.start_time = start_time

        x = np.arange(lon[0], lon[1], dxdy)
        y = np.arange(lat[0], lat[1], dxdy)
        self.lons, self.lats = np.meshgrid(x,y)

        self.path = f.create_directory(date)
        self.thredds_files = f.files_from_thredds(date)
        self.number_of_sims = int(hours/(self.duration/(60*60)))
    
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
        for i in range(24):
            self.inititate_advection_one_member(i)

    def inititate_advection_one_member(self, member):
        #Run for a single member of the current day
        corr = f.correct_file(self.thredds_files, member)
        file = corr[0]
        _member = corr[1]
        
        values = self.set_up_dictionary(self.number_of_sims)
        for i in range(self.number_of_sims):
            for j in range(2):
                o = OceanDrift(loglevel=20)
                r = reader_netCDF_CF_generic.Reader(file, ensemble_member=_member)
                o.add_reader(r)
                if j == 0:
                    o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(hours=i*self.duration)+timedelta(hours=self.start_time))
                    o.run(duration=timedelta(seconds=self.duration), time_step=self.time_step, time_step_output=self.duration, outfile=f'{self.path}/m{member}_f.nc')
                if j == 1:
                    o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(hours=i*self.duration)+timedelta(hours=self.start_time))
                    o.run(duration=timedelta(hours=self.duration), time_step=-self.time_step, time_step_output=-self.duration, outfile=f'{self.path}/m{member}_b.nc')
            _f = xr.open_dataset(f'{self.path}/m{member}_f.nc')
            _b = xr.open_dataset(f'{self.path}/m{member}_b.nc')
            os.remove(f'{self.path}/m{member}_f.nc')
            os.remove(f'{self.path}/m{member}_b.nc')
            values['time'][i+1]['f']['lon'] = np.array(_f['lon'][:,-1])
            values['time'][i+1]['f']['lat'] = np.array(_f['lat'][:,-1])
            values['time'][i+1]['b']['lon'] = np.array(_b['lon'][:,-1])
            values['time'][i+1]['b']['lat'] = np.array(_b['lat'][:,-1])
        np.save(f'{self.path}/values_m{member}.npy', values)
        
        return f'{self.path}/values_m{member}.npy'

    def aggregated(self, start_time, hours):
        file = 'https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_zdepth_be'
        start_time = datetime(start_time[0], start_time[1], start_time[2], start_time[3])

        duration = timedelta(hours=hours)
        values = {'duration': duration, 'time_step': self.time_step, 'dxdy': self.dxdy, 'initial_lon': self.lons, 'initial_lat': self.lats, 'time': {}}
        values['time'][1] = {}
        values['time'][1]['f'] = {}
        values['time'][1]['b'] = {}

        for j in range(2):
            o = OceanDrift(loglevel=20)
            r = reader_netCDF_CF_generic.Reader(file)
            o.add_reader(r)
            if j == 0:
                o.seed_elements(self.lons.ravel(), self.lats.ravel(), time = start_time)
                o.run(duration=duration, time_step=self.time_step, time_step_output=duration, outfile='f.nc')
            if j == 1:
                o.seed_elements(self.lons.ravel(), self.lats.ravel(), time = start_time)
                o.run(duration=duration, time_step=-self.time_step, time_step_output=duration, outfile='b.nc')
        _f = xr.open_dataset('f.nc')
        _b = xr.open_dataset('b.nc')
        os.remove('f.nc')
        os.remove('b.nc')
        values['time'][1]['f']['lon'] = np.array(_f['lon'][:,-1])
        values['time'][1]['f']['lat'] = np.array(_f['lat'][:,-1])
        values['time'][1]['b']['lon'] = np.array(_b['lon'][:,-1])
        values['time'][1]['b']['lat'] = np.array(_b['lat'][:,-1])
        np.save('values.npy', values)

    def aggregated_backwards(self, start_time, at_time = 0, duration=1, par=False, outfile='values'):
        file = 'https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_zdepth_be'
        start_time = datetime(start_time[0], start_time[1], start_time[2], start_time[3])

        o = OceanDrift(loglevel=20)
        r = reader_netCDF_CF_generic.Reader(file)
        o.add_reader(r)
        o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=start_time+timedelta(hours = at_time))
        o.run(duration=timedelta(hours=duration), time_step=-self.time_step, time_step_output=timedelta(hours=duration), outfile=f'{outfile}.nc')

        _b = xr.open_dataset(f'{outfile}.nc')
        os.remove(f'{outfile}.nc')
        values = {'duration': duration, 'time_step': self.time_step, 'dxdy': self.dxdy, 'initial_lon': self.lons, 'initial_lat': self.lats, 
        'time': {1:{'f': {}, 'b':{}}}}
        values['time'][1]['b']['lon'] = np.array(_b['lon'][:,-1])
        values['time'][1]['b']['lat'] = np.array(_b['lat'][:,-1])
        values['time'][1]['f']['lon'] = np.zeros_like(_b['lon'][:,-1])
        values['time'][1]['f']['lat'] = np.zeros_like(_b['lat'][:,-1])

        np.save('values.npy', values)
        if par == True:
            path = f'tmp/{outfile}.npy'
            np.save(path, values)  
        else:
            path = f'{outfile}.npy'
            np.save(path)
        return path
    
if __name__ == '__main__':

    r = particle_grid_displacement([12.5, 16], [68,70], 900, 0.01, 2, '20220529', 12)
    r.displacement_one_member(1)

    """
    lon = [8.1,8.4]
    lat = [70,70.3]

    time_step = 900
    duration = 7200
    dxdy = 0.1
    start_time = 12
    #initiate_file_all_members(lon, lat, time_step, dxdy, duration, '20220420')
    r = initiate_particle_advection(lon,lat,time_step,dxdy,duration,'20220529', hours=2, start_time=start_time)
    for i in range(6):
        r.inititate_advection_one_member(i)
    #r.initiate_advection_all_members
    #r.inititate_advection_one_member(0)
    #r.aggregated([2022,5,30,12], hours = 2)
    """
#https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_eps_files/barents_eps_20220420T18Z.nc
#https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_eps_files/barents_eps_20220420T12Z.nc
#https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_eps_files/barents_eps_20220420T06Z.nc
#https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_eps_files/barents_eps_20220420T00Z.nc