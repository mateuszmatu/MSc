#!/usr/bin/env python
from datetime import datetime, timedelta

from opendrift.readers import reader_double_gyre
from opendrift.models.oceandrift import OceanDrift
import numpy as np
import os
import xarray as xr
import re
def example_double_gyre(at_time = 3, duration = 15, time_step = 0.5, dxdy = 0.02, epsilon=0.25, omega=0.682, A=0.1):
    if not os.path.exists('flow_maps'):
        os.mkdir('flow_maps')
    if not os.path.exists('flow_maps/double_gyre'):
        os.mkdir('flow_maps/double_gyre')
    save_str = f'{at_time}_{duration}_{time_step}_{dxdy}_{epsilon}_{omega}_{A}'
    save_str = re.sub('\.', '', save_str)
    path_to_file = f'flow_maps/double_gyre/{save_str}.npy'    
    for i in range(2):
        o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('drift:advection_scheme', 'runge-kutta4')

        double_gyre = reader_double_gyre.Reader(epsilon=epsilon, omega=omega, A=A)
        print (double_gyre)

        o.add_reader(double_gyre)

        proj = double_gyre.proj
        x = np.arange(double_gyre.xmin, double_gyre.xmax, dxdy)
        y = np.arange(double_gyre.ymin, double_gyre.ymax, dxdy)

        X, Y = np.meshgrid(x,y)
        lons, lats = proj(X,Y, inverse=True)
        if i == 0:
            o.seed_elements(lons.ravel(), lats.ravel(),
                            time=double_gyre.initial_time+timedelta(seconds=at_time))
            o.run(duration=timedelta(seconds=duration), time_step=time_step, time_step_output=duration, outfile='flow_maps/double_gyre/double_gyre_f.nc')
        elif i == 1:
            o.seed_elements(lons.ravel(), lats.ravel(),
                            time=double_gyre.initial_time+timedelta(seconds=duration+at_time))
            o.run(duration=timedelta(seconds=duration), time_step=-time_step, time_step_output=duration, outfile='flow_maps/double_gyre/double_gyre_b.nc')

    _f = xr.open_dataset('flow_maps/double_gyre/double_gyre_f.nc')
    _b = xr.open_dataset('flow_maps/double_gyre/double_gyre_b.nc')
    os.remove('flow_maps/double_gyre/double_gyre_f.nc')
    os.remove('flow_maps/double_gyre/double_gyre_b.nc')
    values = {'duration': duration, 'time_step': time_step, 'dxdy': dxdy, 'initial_lon': lons, 'initial_lat': lats, 
        'time': {1:{'f': {}, 'b':{}}}}
    values['time'][1]['f']['lon'] = np.array(_f['lon'][:,-1])
    values['time'][1]['f']['lat'] = np.array(_f['lat'][:,-1])
    values['time'][1]['b']['lon'] = np.array(_b['lon'][:,-1])
    values['time'][1]['b']['lat'] = np.array(_b['lat'][:,-1])
    np.save(f'flow_maps/double_gyre/{save_str}.npy', values)

    return path_to_file

def example_double_gyre_backward(at_time = 3, duration = 15, time_step = 0.5, dxdy = 0.02, epsilon=0.25, omega=0.682, A=0.1):
    if not os.path.exists('flow_maps'):
        os.mkdir('flow_maps')
    if not os.path.exists('flow_maps/double_gyre'):
        os.mkdir('flow_maps/double_gyre')
    save_str = f'{at_time}_{duration}_{time_step}_{dxdy}_{epsilon}_{omega}_{A}'
    save_str = re.sub('\.', '', save_str)
    path_to_file = f'flow_maps/double_gyre/b_{save_str}.npy' 

    o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
    o.set_config('environment:fallback:land_binary_mask', 0)
    o.set_config('drift:advection_scheme', 'runge-kutta4')

    double_gyre = reader_double_gyre.Reader(epsilon=epsilon, omega=omega, A=A)
    print (double_gyre)

    o.add_reader(double_gyre)

    proj = double_gyre.proj
    x = np.arange(double_gyre.xmin, double_gyre.xmax, dxdy)
    y = np.arange(double_gyre.ymin, double_gyre.ymax, dxdy)

    X, Y = np.meshgrid(x,y)
    lons, lats = proj(X,Y, inverse=True)
    o.seed_elements(lons.ravel(), lats.ravel(),
                    time=double_gyre.initial_time+timedelta(seconds=duration+at_time))
    o.run(duration=timedelta(seconds=duration), time_step=-time_step, time_step_output=duration, outfile='flow_maps/double_gyre/double_gyre_b.nc')

    _b = xr.open_dataset('flow_maps/double_gyre/double_gyre_b.nc')
    os.remove('flow_maps/double_gyre/double_gyre_b.nc')
    values = {'duration': duration, 'time_step': time_step, 'dxdy': dxdy, 'initial_lon': lons, 'initial_lat': lats, 
        'time': {1:{'f': {}, 'b':{}}}}
    values['time'][1]['b']['lon'] = np.array(_b['lon'][:,-1])
    values['time'][1]['b']['lat'] = np.array(_b['lat'][:,-1])
    values['time'][1]['f']['lon'] = np.zeros_like(_b['lon'][:,-1])
    values['time'][1]['f']['lat'] = np.zeros_like(_b['lat'][:,-1])

    np.save(f'flow_maps/double_gyre/b_{save_str}.npy', values)

    return path_to_file
if __name__ == '__main__':
    example_double_gyre()
    