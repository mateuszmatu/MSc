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
from joblib import Parallel, delayed

class particle_grid_displacement:
    def __init__(self, lon, lat, time_step, sep, dur, date=None, hours=None):
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

        if date is not None:
            self.path = f.create_directory(date)
            self.tf = f.files_from_thredds(date, 'thredds_urls.txt')
            self.obf = f.files_from_thredds(date, 'old_barents.txt')[0]
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

        hour_adjustment = f.hour_adjustment(member)

        if self.h is None:
            hours = hour_adjustment
        elif self.h is not None:
            hours = hour_adjustment + self.h

        corr = f.correct_file(self.tf, member)
        file = corr[0]
        _member = corr[1]

        values = self._set_up_dictionary
        if outfile is None:
            string = f'{self.path}/{self.name}_h{self.h}_m{member}'
        elif outfile is not None:
            string = f'{outfile}'
        
        for i in range(2):
            o = OceanDrift(loglevel=30)
            r = reader_netCDF_CF_generic.Reader(file, ensemble_member=_member)
            o.add_reader(r)
            if i == 0 and forwards is True:
                o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(hours=hours))
                o.run(duration=timedelta(hours=self.dur), time_step=timedelta(seconds=self.ts), time_step_output=timedelta(hours=self.dur), outfile=f'{string}_f.nc')
            if i == 1 and backwards is True:
                o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(hours=hours))
                o.run(duration=timedelta(hours=self.dur), time_step=timedelta(seconds=-self.ts), time_step_output=timedelta(hours=-self.dur), outfile=f'{string}_b.nc')

        if forwards is True:
            _f = xr.open_dataset(f'{string}_f.nc')
            os.remove(f'{string}_f.nc')
            values['f']['lon'] = np.array(_f['lon'][:,-1])
            values['f']['lat'] = np.array(_f['lat'][:,-1])
            if backwards is False:
                values['b']['lon'] = np.zeros_like(_f['lon'][:,-1])
                values['b']['lat'] = np.zeros_like(_f['lat'][:,-1])

        if backwards is True:
            _b = xr.open_dataset(f'{string}_b.nc')
            os.remove(f'{string}_b.nc')
            values['b']['lon'] = np.array(_b['lon'][:,-1])
            values['b']['lat'] = np.array(_b['lat'][:,-1])
            if forwards is False:
                values['f']['lon'] = np.zeros_like(_b['lon'][:,-1])
                values['f']['lat'] = np.zeros_like(_b['lat'][:,-1])
        
        string = string + '.npy'
        np.save(string, values)

        return string
        
    def displacement_all_members(self, num_of_members = 24, forwards = True, backwards = True, njobs = 2):
        '''
            Runs the displacement_one_member() function for all members in the eps.
        Args:
            num_of_members   [int]  : number of members in the eps. Default 24 ---> number of members in the barents 2.5 eps.
            forwards         [bool] : True if particles should be advected forwards in time, False otherwise
            backwards        [bool] : True if particles should be advected backwards in time, False otherwise
            njobs            [int]  : How many jobs should run at the same time. Default = 2
        Returns:
            files            [list] : A list of filepath to each of the member files   
        '''
        """
        for i in range(24):
            addhour = f.hour_adjustment(i)
            file = self.displacement_one_member(i, forwards = forwards, backwards = backwards, hour_adjustment = addhour)
        """
        def par(i):
            file = self.displacement_one_member(i, forwards = forwards, backwards = backwards)
            return file
        files = Parallel(n_jobs=njobs)(delayed(par)(i) for i in range(num_of_members))
        return files
    
    def aggregated_displacement(self, start_time, at_time, outfile, par=True):
        '''
            Seeds particles using an aggregated dataset, so that they can be advected for a longer time.
        Args:   
            start_time  [list]  :   a list of the date [year, month, day, hour]
            outfile     [str]   :   name of outfile
            par         [bool]  :   if this is used in a method that uses a parallel loop (not really important)
        Returns:
            path        [str]   :   a string with path to file
        '''

        file = 'https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_zdepth_be'
        st = datetime(start_time[0], start_time[1], start_time[2], start_time[3])

        o = OceanDrift(loglevel=30)
        r = reader_netCDF_CF_generic.Reader(file)
        o.add_reader(r)
        o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=st+timedelta(hours=at_time))
        o.run(duration=timedelta(hours=self.dur), time_step=-self.ts, time_step_output=timedelta(hours=self.dur), outfile=f'{outfile}.nc')

        _b = xr.open_dataset(f'{outfile}.nc')
        os.remove(f'{outfile}.nc')

        values = self._set_up_dictionary

        values['b']['lon'] = np.array(_b['lon'][:,-1])
        values['b']['lat'] = np.array(_b['lat'][:,-1])
        values['f']['lon'] = np.zeros_like(_b['lon'][:,-1])
        values['f']['lat'] = np.zeros_like(_b['lat'][:,-1])

        if par is True:
            path = f'tmp/{outfile}.npy'
            np.save(path, values)
        elif par is False:
            path = f'{outfile}.npy'
            np.save(path)
        return path
    
    def non_eps_advection(self, dur, outfile):
        
        o = OceanDrift(loglevel=30)
        r = reader_netCDF_CF_generic.Reader(self.obf)
        o.add_reader(r)
        o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(hours=self.h)-timedelta(hours=dur))
        o.run(duration=timedelta(hours=dur), time_step=timedelta(seconds=self.ts), time_step_output=timedelta(hours=dur), outfile=outfile)

        return outfile

 
    
if __name__ == '__main__':

    r = particle_grid_displacement([12.5, 16], [68,70], 900, 0.01, 2, '20220829')
    r.displacement_all_members()

