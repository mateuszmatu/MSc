#!/usr/bin/env python
from datetime import datetime, timedelta

from opendrift.readers import reader_double_gyre
from opendrift.models.oceandrift import OceanDrift
import numpy as np
import os
import xarray as xr
import re
import matplotlib.pyplot as plt

def example_double_gyre(at_time = 3, dur = 15, time_step = 0.5, sep = 0.01, epsilon = 0.25, omega = 0.682, A = 0.1, forwards = True, backwards = True, outfile = None):
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

        if i == 0 and forwards is True:
            o.seed_elements(lons.ravel(), lats.ravel(),
                            time=double_gyre.initial_time+timedelta(seconds=at_time))
            o.run(duration=timedelta(seconds=dur), time_step=time_step, time_step_output=dur, outfile=f'flow_maps/double_gyre/{outfile}_f.nc')
        elif i == 1 and backwards is True:
            o.seed_elements(lons.ravel(), lats.ravel(),
                            time=double_gyre.initial_time+timedelta(seconds=at_time))
            o.run(duration=timedelta(seconds=dur), time_step=-time_step, time_step_output=dur, outfile=f'flow_maps/double_gyre/{outfile}_b.nc')
    
    values = {'duration': dur, 'time_step': time_step, 'in_separation': sep, 'in_lons': lons, 'in_lats': lats, 'f': {}, 'b': {}}

    if forwards is True:
        _f = xr.open_dataset(f'flow_maps/double_gyre/{outfile}_f.nc')
        os.remove(f'flow_maps/double_gyre/{outfile}_f.nc')
        values['f']['lon'] = np.array(_f['lon'][:,-1])
        values['f']['lat'] = np.array(_f['lat'][:,-1])
        if backwards is False:
            values['b']['lon'] = np.zeros_like(_f['lon'][:,-1])
            values['b']['lat'] = np.zeros_like(_f['lat'][:,-1])

    if backwards is True:
        _b = xr.open_dataset(f'flow_maps/double_gyre/{outfile}_b.nc')
        os.remove(f'flow_maps/double_gyre/{outfile}_b.nc')
        values['b']['lon'] = np.array(_b['lon'][:,-1])
        values['b']['lat'] = np.array(_b['lat'][:,-1])
        if forwards is False:
            values['f']['lon'] = np.zeros_like(_b['lon'][:,-1])
            values['f']['lat'] = np.zeros_like(_b['lat'][:,-1])
    
    np.save(f'flow_maps/double_gyre/{save_str}.npy', values)
    
    return path

def double_gyre_velocity_field(at_time = 3, epsilon = 0.25, omega = 0.682, A = 0.1):
    '''
        A function for showing the velocity field in a double gyre.
    Args:
        at_time     [int]   :   the time (seconds) at which particles should be released
        epsilon     [float] :   no clue
        omega       [float] :   no clue
        A           [float] :   no clue
    Returns:
        x           [arr]   :   an array of uniformly spaced values from 0 to 2
        y           [arr]   :   an array of uniformly spaced values from 0 to 1
        u           [arr]   :   an array containing the u velocities
        v           [arr]   :   an array containing the v velocities
        infor_str   [list]  :   a list containing the parameters used
    '''
    a = epsilon*np.sin(omega*at_time)
    b = 1-2*epsilon*np.sin(omega*at_time)

    def f(x):
        return a*x*x+b*x
    
    def fx(x):
        return 2*x*a+b
    
    x,y = np.meshgrid(np.linspace(0,2,20), np.linspace(0,1,20))
    
    u = -np.pi*A*np.sin(np.pi*f(x))*np.cos(np.pi*y)
    v = np.pi*A*np.cos(np.pi*f(x))*np.sin(np.pi*y)*fx(x)

    info_str = [at_time, epsilon, omega, A]
    return x,y,u,v,info_str

if __name__ == '__main__':
    from LCS import FTLE
    LCS1 = FTLE(example_double_gyre(epsilon=0.25, A=0.1, omega=2*np.pi/10, at_time=0.5, forwards=False))
    LCS2 = FTLE(example_double_gyre(epsilon=0.25, A=0.1, omega=0.7, at_time=0.5, forwards=False))
    LCS3 = FTLE(example_double_gyre(epsilon=0.25, A=0.2, omega=2*np.pi/10, at_time=0.5, forwards=False))

    ex_1 = double_gyre_velocity_field(epsilon=0.25, A=0.1, omega=2*np.pi/10, at_time=0.5)
    ex_2 = double_gyre_velocity_field(epsilon=0.25, A=0.1, omega=0.7, at_time=0.5)
    ex_3 = double_gyre_velocity_field(epsilon=0.25, A=0.2, omega=2*np.pi/10, at_time=0.5)

    fig, ax = plt.subplots(2,3, figsize=(6,12))
    #ax[0].pcolormesh(LCS['lon'], LCS['lat'], LCS['ALCS'], cmap='jet')
    Q1 = ax[0,0].quiver(ex_1[0], ex_1[1], ex_1[2], ex_1[3], scale=10)
    Q2 = ax[0,1].quiver(ex_2[0], ex_2[1], ex_2[2], ex_2[3], scale=10)
    Q3 = ax[0,2].quiver(ex_3[0], ex_3[1], ex_3[2], ex_3[3], scale=10)

    ax[1,0].pcolormesh(LCS1['lon'], LCS1['lat'], LCS1['ALCS'], cmap='jet')
    ax[1,1].pcolormesh(LCS2['lon'], LCS2['lat'], LCS2['ALCS'], cmap='jet')
    ax[1,2].pcolormesh(LCS3['lon'], LCS3['lat'], LCS3['ALCS'], cmap='jet')

    qk1 = ax[0,0].quiverkey(Q1, 1.05, 1.05, 1, r'$1 \frac{m}{s}$', labelpos='E', coordinates='axes')
    qk2 = ax[0,1].quiverkey(Q2, 1.05, 1.05, 1, r'$1 \frac{m}{s}$', labelpos='E', coordinates='axes')
    qk3 = ax[0,2].quiverkey(Q3, 1.05, 1.05, 1, r'$1 \frac{m}{s}$', labelpos='E', coordinates='axes')

    ax[0,0].set_title(f'time: {ex_1[4][0]}, epsilon: {ex_1[4][1]}, omega: {ex_1[4][2]:.3f}, A: {ex_1[4][3]}')
    ax[0,1].set_title(f'time: {ex_2[4][0]}, epsilon: {ex_2[4][1]}, omega: {ex_2[4][2]:.3f}, A: {ex_2[4][3]}')
    ax[0,2].set_title(f'time: {ex_3[4][0]}, epsilon: {ex_3[4][1]}, omega: {ex_3[4][2]:.3f}, A: {ex_3[4][3]}')
    plt.show()