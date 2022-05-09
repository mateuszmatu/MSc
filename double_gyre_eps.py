from example_double_gyre import example_double_gyre as edg
from example_double_gyre import example_double_gyre_backward as edgb
from LCS import FTLE
from opendrift.readers import reader_double_gyre
from opendrift.models.oceandrift import OceanDrift
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from datetime import datetime, timedelta
import xarray as xr

def ALCS_plot(LCS):
    plt.pcolormesh(LCS['lon'], LCS['lat'], LCS['ALCS'][1], cmap='jet')
    plt.title('ALCS')
    plt.colorbar()
    plt.show()

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

def animation(i):
    file = edgb(times[i],15,0.5,0.02,0.25,0.682,0.1)
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

    double_gyre = reader_double_gyre.Reader(epsilon=.1, omega=0.628, A=0.1)
    print (double_gyre)

    o.add_reader(double_gyre)

    x = [.9]
    y = [.5]
    lon, lat = double_gyre.xy2lonlat(x, y)
    time_step = times[1]-times[0]
    o.seed_elements(lon, lat, radius=.1, number=5000,
                    time=double_gyre.initial_time+timedelta(seconds=18))

    o.run(duration=timedelta(seconds=20), time_step=time_step, outfile='flow_maps/double_gyre/particle_advection.nc')
    
def LCS_and_particles_advection_animation(outfile, frames, t0, t1, duration = 15, time_step = 0.5, dxdy = 0.02, epsilon=0.25, omega=0.682, A=0.1):

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
                        time=double_gyre.initial_time+timedelta(seconds=18))

        o.run(duration=timedelta(seconds=20), time_step=time_step, outfile='flow_maps/double_gyre/particle_advection.nc')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    times = np.linspace(t0,t1,frames)
    advect_particles()
    pos = xr.load_dataset('flow_maps/double_gyre/particle_advection.nc')

    animation = FuncAnimation(fig, func=animation, interval=200, frames=frames)
    animation.save(f'gifs/{outfile}.gif', dpi = 80, writer='imagemagick')

LCS_and_particles_advection_animation('LCS3', 50, 1, 20, omega=0.1)