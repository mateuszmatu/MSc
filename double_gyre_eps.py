
from example_double_gyre import example_double_gyre_backward as edgb
from example_double_gyre import double_gyre_for_animation as dgfa
from LCS import FTLE
from opendrift.readers import reader_double_gyre
from opendrift.models.oceandrift import OceanDrift
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from datetime import datetime, timedelta
import xarray as xr
import time
import os
from joblib import Parallel, delayed


def ALCS_plot(LCS, i):
    plt.pcolormesh(LCS['lon'], LCS['lat'], LCS['ALCS'][1], cmap='jet')
    plt.title('ALCS')
    plt.colorbar()
    plt.savefig(f'double_gyre_example_{i}.png')
    plt.clf()

def difference(LCS, LCS2):
    fig, ax = plt.subplots(3)
    a = ax[0].pcolormesh(LCS['lon'], LCS['lat'], LCS['ALCS'][1], cmap='jet')
    b = ax[1].pcolormesh(LCS2['lon'], LCS2['lat'], LCS2['ALCS'][1], cmap='jet')
    diff = LCS['ALCS'][1]-LCS2['ALCS'][1]
    c = ax[2].pcolormesh(LCS['lon'], LCS['lat'], diff, cmap='jet')

    ax[0].set(title='ALCS 1')
    ax[1].set(title='ALCS 2')
    ax[2].set(title='ALCS 1 - ALCS 2')
    l = [a,b,c]
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    for i in range(3):
        divider = make_axes_locatable(ax[i])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(l[i], cax=cax)
    plt.show()

def LCS_times_series(time):
    time_series = {}
    for t in time:
        file = edgb(t,15,0.5,0.02,0.1,0.682,0.1)
        LCS = FTLE(file)
        time_series[t] = LCS['ALCS'][1]
    
    return time_series
    
def LCS_and_particles_advection_animation_gyre(outfile, frames, t0, t1, duration = 15, time_step = 0.5, dxdy = 0.02, epsilon=0.25, omega=0.682, A=0.1):

    def animation(i):
        file = edgb(times[i],duration,time_step,dxdy,epsilon,omega,A)
        LCS = FTLE(file)
        x = []
        x = LCS['ALCS'][1]
        posx = []
        posy = []
        posx = pos['lon'][:,i]
        posy = pos['lat'][:,i]
        ax.clear()
        ax.set_title(f'frame {i}/{frames}')
        ax.pcolormesh(LCS['lon'], LCS['lat'], x, cmap='jet')
        ax.scatter(posx, posy, s=5, color='black')

    def advect_particles():
        o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('drift:advection_scheme', 'runge-kutta4')

        double_gyre = reader_double_gyre.Reader(epsilon=epsilon, omega=omega, A=A)
        print (double_gyre)

        o.add_reader(double_gyre)

        x = [.9]
        y = [.5]
        lon, lat = double_gyre.xy2lonlat(x, y)
        time_step = times[1]-times[0]
        o.seed_elements(lon, lat, radius=.1, number=5000,
                        time=double_gyre.initial_time+timedelta(seconds=t0))
        o.run(duration=timedelta(seconds=40), time_step=time_step, outfile='flow_maps/double_gyre/particle_advection.nc')

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    times = np.linspace(t0,t1,frames)
    #advect_particles()
    pos = xr.load_dataset('flow_maps/double_gyre/particle_advection.nc')

    ani = Parallel(n_jobs=8)(delayed(FuncAnimation(fig, func=animation, interval=200, frames=frames)))
    #animation = FuncAnimation(fig, func=animation, interval=200, frames=frames)
    ani.save(f'gifs/{outfile}.gif', dpi = 80, writer='imagemagick')

def LCS_and_particles_advection_animation_gyre_par(outfile, frames, t0, t1, duration = 15, time_step = 0.5, dxdy = 0.02, epsilon=0.25, omega=0.682, A=0.1):

    def animation(i):
        LCS = arr[i]
        x = []
        x = LCS
        posx = []
        posy = []
        posx = pos['lon'][:,i]
        posy = pos['lat'][:,i]
        ax.clear()
        ax.set_title(f'frame {i}/{frames}')
        ax.pcolormesh(in_x, in_y, x, cmap='jet')
        ax.scatter(posx, posy, s=5, color='black')

    def advect_particles():
        o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('drift:advection_scheme', 'runge-kutta4')

        double_gyre = reader_double_gyre.Reader(epsilon=epsilon, omega=omega, A=A)
        print (double_gyre)

        o.add_reader(double_gyre)

        x = [.9]
        y = [.5]
        lon, lat = double_gyre.xy2lonlat(x, y)
        time_step = times[1]-times[0]
        o.seed_elements(lon, lat, radius=.1, number=5000,
                        time=double_gyre.initial_time+timedelta(seconds=t0))
        o.run(duration=timedelta(seconds=t1), time_step=time_step, outfile='flow_maps/double_gyre/particle_advection.nc')

    def par(i):
        file = dgfa(times[i],duration,time_step,dxdy,epsilon,omega,A,outfile=i)
        LCS = FTLE(file)
        return LCS['ALCS'][1]
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    times = np.linspace(t0,t1,frames)
    advect_particles()
    pos = xr.load_dataset('flow_maps/double_gyre/particle_advection.nc')

    arr = Parallel(n_jobs=4)(delayed(par)(i) for i in range(frames))
    data = np.load('flow_maps/double_gyre/1.npy',allow_pickle='TRUE').item()
    in_x = data['initial_lon']
    in_y = data['initial_lat']
    ani = FuncAnimation(fig, func=animation, interval=200, frames=frames)
    ani.save(f'gifs/{outfile}.gif', dpi = 80, writer = 'imagemagick')

    def remove_files(i):
        os.remove(f'flow_maps/double_gyre/{i}.npy')

    Parallel(n_jobs=4)(delayed(remove_files)(i) for i in range(frames))


    
"""
file = edgb(at_time = 3, duration = 15, time_step = 0.5, dxdy = 0.02, epsilon=0.25, omega=0.682, A=0.1)
LCS1 = FTLE(file)
file = edgb(at_time = 3, duration = 15, time_step = 0.5, dxdy = 0.02, epsilon=0.20, omega=0.650, A=0.05)
LCS2 = FTLE(file)
file = edgb(at_time = 3, duration = 15, time_step = 0.5, dxdy = 0.02, epsilon=0.30, omega=0.700, A=0.15)
LCS3 = FTLE(file)
LCS_arr = [LCS1, LCS2, LCS3]
for i, j in zip(LCS_arr, range(3)):
    ALCS_plot(i,j)
"""

start = time.perf_counter()
LCS_and_particles_advection_animation_gyre_par('para', 100, 1, 40)
end = time.perf_counter()
print(end-start)