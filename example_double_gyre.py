#!/usr/bin/env python
from datetime import datetime, timedelta

from opendrift.readers import reader_double_gyre
from opendrift.models.oceandrift import OceanDrift
import numpy as np
import os
import xarray as xr
import re

def example_double_gyre(at_time = 3, dur = 15, time_step = 0.5, sep = 0.02, epsilon = 0.25, omega = 0.682, A = 0.1, forwards = True, backwards = True, outfile = None):
    '''
        A script for setting up a controlled double gyre environment, and then advecting particles in this environment to calculate LCS later.
    Args:
        at_time     [int]   :   the time (seconds) at which particles should be released
        dur         [int]   :   the duration (seconds) of the simulation
        time_step   [float] :   the time step (seconds) of the simulation
        sep         [float] :   the initial separation between seeded particles
        epsilon     [float] :   no clue
        omega       [float] :   no clue
        A           [float] :   no clue
        forwards    [bool]  :   True if particles should be advected forwards in time, False otherwise
        backwards   [bool]  :   True if particles should be advected backwards in time, False otherwise
        outfile     [str]   :   the name of the saved file
    Returns:
        path        [str]   :   a path to the saved file containing results of simulation
    '''

    if not os.path.exists('flow_maps'):
        os.mkdir('flow_maps')
    if not os.path.exists('flow_maps/double_gyre'):
        os.mkdir('flow_maps/double_gyre')
    
    if outfile is None:
        save_str = f'{at_time}_{dur}_{time_step}_{sep}_{epsilon}_{omega}_{A}'
        save_str = re.sub('\.', '-',save_str)
    elif outfile is not None:
        save_str = outfile
    path = f'flow_maps/double_gyre/{save_str}.npy'

    for i in range(2):
        o = OceanDrift(loglevel=20)
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('drift:advection_scheme', 'runge-kutta4')
        double_gyre = reader_double_gyre.Reader(epsilon=epsilon, omega=omega, A=A)
        o.add_reader(double_gyre)

        proj = double_gyre.proj
        x = np.arange(double_gyre.xmin, double_gyre.xmax, sep)
        y = np.arange(double_gyre.ymin, double_gyre.ymax, sep)

        X, Y = np.meshgrid(x,y)
        lons, lats = proj(X,Y, inverse=True)

        if i == 0:
            o.seed_elements(lons.ravel(), lats.ravel(),
                            time=double_gyre.initial_time+timedelta(seconds=at_time))
            o.run(duration=timedelta(seconds=dur), time_step=time_step, time_step_output=dur, outfile='flow_maps/double_gyre/double_gyre_f.nc')
        elif i == 1:
            o.seed_elements(lons.ravel(), lats.ravel(),
                            time=double_gyre.initial_time+timedelta(seconds=at_time))
            o.run(duration=timedelta(seconds=dur), time_step=-time_step, time_step_output=dur, outfile='flow_maps/double_gyre/double_gyre_b.nc')
    
    values = {'duration': dur, 'time_step': time_step, 'in_separation': sep, 'in_lons': lons, 'in_lats': lats, 'f': {}, 'b': {}}

    if forwards is True:
        _f = xr.open_dataset('flow_maps/double_gyre/double_gyre_f.nc')
        os.remove('flow_maps/double_gyre/double_gyre_f.nc')
        values['f']['lon'] = np.array(_f['lon'][:,-1])
        values['f']['lat'] = np.array(_f['lat'][:,-1])
        if backwards is False:
            values['b']['lon'] = np.zeros_like(_f['lon'][:,-1])
            values['b']['lat'] = np.zeros_like(_f['lat'][:,-1])

    if backwards is True:
        _b = xr.open_dataset('flow_maps/double_gyre/double_gyre_b.nc')
        os.remove('flow_maps/double_gyre/double_gyre_b.nc')
        values['b']['lon'] = np.array(_b['lon'][:,-1])
        values['b']['lat'] = np.array(_b['lat'][:,-1])
        if forwards is False:
            values['f']['lon'] = np.zeros_like(_b['lon'][:,-1])
            values['f']['lat'] = np.zeros_like(_b['lat'][:,-1])
    
    np.save(f'flow_maps/double_gyre/{save_str}.npy', values)

    return path

if __name__ == '__main__':
    example_double_gyre()
    