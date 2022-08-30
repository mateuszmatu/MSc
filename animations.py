from example_double_gyre import example_double_gyre as edg
from grid_advection import particle_grid_displacement as pgd
from LCS import FTLE
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import cartopy
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1 import make_axes_locatable
from joblib import Parallel, delayed
from datetime import datetime, timedelta
import time
import os
from opendrift.readers import reader_double_gyre
from opendrift.models.oceandrift import OceanDrift

class animations:
    def __init__(self, lon=None, lat=None, time_step=None, sep=None, dur=None, date=None, hours=None, frames=5, jobs=2):
        '''
            Class for creating different animations of LCS and particles
        Args:
            lon         [list]  :   list of min and max lon
            lat         [list]  :   list of min and max lat
            time_step   [int]   :   time step in seconds
            sep         [float] :   initial separation between particles 
            dur         [int]   :   duration for how long particles should drift
            date        [str]   :   the date for the dataset
            hours       [int]   :   hours after dataset start
            frames      [int]   :   amount of frames in animation
            jobs        [int]   :   amount of loops running at the same time in the parallel loop
        '''

        self.lon = lon
        self.lat = lat
        self.ts = time_step
        self.sep = sep
        self.dur = dur
        #self.date = date       these two are not used as of now
        #self.h = hours
        self.frames = frames
        self.jobs = jobs
    
    def LCS_and_particle_advection(self, start_time, outfile, LCS_start_time=0):
        '''
            Creates an animation of particle advection and LCS field.
        Args:
            start_time      [list]  :   a list of the date [year, month, day, hour]
            outfile         [str]   :   name of the saved gif
            LCS_start_time  [int]   :   start time of LCS field compared to particle advection simulation (hours)
        '''

        def par(i):
            ip = pgd(self.lon, self.lat, self.ts, self.sep, self.dur)
            file = ip.aggregated_displacement(start_time, at_time=i, outfile=i, par=True)
            LCS = FTLE(file)
            return LCS['ALCS']
        
        def animation(i):
            x = []
            x = arr[i]
            posx = []
            posy = []
            posx = pos['lon'][:,i+LCS_start_time]
            posy = pos['lat'][:,i+LCS_start_time]
            ax.set_title(f'frame {i}/{self.frames}')
            ax.pcolormesh(in_x, in_y, x, cmap='jet', vmin=-0.2, vmax=0.5, transform=plate)
            ax.scatter(posx, posy, s=5, color='black', transform=plate)
            ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
            ax.coastlines()
            ax.gridlines(draw_labels=True)

        def remove_files(i):
            os.remove(f'tmp/{i}.npy')
        
        transform = ccrs.NorthPolarStereo()
        plate = ccrs.PlateCarree()
        fig = plt.figure(figsize=(20,20))
        ax = plt.axes(projection=transform)
        pos = xr.load_dataset('particle_drift.nc')

        arr = Parallel(n_jobs=self.jobs)(delayed(par)(i) for i in range(self.frames))
        data = np.load('tmp/0.npy', allow_pickle='True').item()
        in_x = data['in_lons']
        in_y = data['in_lats']
        ani=FuncAnimation(fig, func=animation, interval=300, frames=self.frames)
        ani.save(f'gifs/{outfile}.gif', dpi=80, writer='imagemagick')

        Parallel(n_jobs=4)(delayed(remove_files)(i) for i in range(self.frames))    #cleanup

    def LCS_and_particle_advection_gyre(self, t0, t1, outfile, epsilon=0.25, omega=0.682, A=0.1):
        '''
            Creates an animation of particle advection and LCS field in the double gyre example
        Args:
            t0          [int]   :   start of simulation (seconds)
            t1          [int]   :   end of simulation (seconds)
            outfile     [str]   :   outfile name
            epsilon     [float] :   No clue
            omega       [float] :   No clue
            A           [float] :   No clue
        '''

        def par(i):
            file = edg(at_time=times[i], dur=self.dur, time_step=self.ts, sep=self.sep, epsilon=epsilon, omega=omega, A=A, forwards=False, outfile=i)
            LCS = FTLE(file)
            return LCS['ALCS']
        
        def advect_particles():
            o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
            o.set_config('environment:fallback:land_binary_mask', 0)
            o.set_config('drift:advection_scheme', 'runge-kutta4')

            double_gyre = reader_double_gyre.Reader(epsilon=epsilon, omega=omega, A=A)

            o.add_reader(double_gyre)

            x = [.9]
            y = [.5]
            lon, lat = double_gyre.xy2lonlat(x, y)
            time_step = times[1]-times[0]
            o.seed_elements(lon, lat, radius=.1, number=5000,
                            time=double_gyre.initial_time+timedelta(seconds=t0))
            o.run(duration=timedelta(seconds=t1), time_step=time_step, outfile='flow_maps/double_gyre/particle_advection.nc')

        def animation(i):
            LCS = arr[i]
            x = []
            x = LCS
            posx = []
            posy = []
            posx = pos['lon'][:,i]
            posy = pos['lat'][:,i]
            ax.clear()
            ax.set_title(f'frame {i}/{self.frames}')
            ax.pcolormesh(in_x, in_y, x, cmap='jet')
            ax.scatter(posx, posy, s=5, color='black')


        def remove_files(i):
            os.remove(f'flow_maps/double_gyre/{i}.npy')

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        times = np.linspace(t0,t1,self.frames)
        advect_particles()
        pos = xr.load_dataset('flow_maps/double_gyre/particle_advection.nc')

        arr = Parallel(n_jobs=self.jobs)(delayed(par)(i) for i in range(self.frames))
        data = np.load('flow_maps/double_gyre/1.npy',allow_pickle='TRUE').item()
        in_x = data['in_lons']
        in_y = data['in_lats']
        ani = FuncAnimation(fig, func=animation, interval=200, frames=self.frames)
        ani.save(f'gifs/{outfile}.gif', dpi = 80, writer = 'imagemagick')

        Parallel(n_jobs=self.jobs)(delayed(remove_files)(i) for i in range(self.frames))

if __name__ == '__main__':
    lon = [14.5,17.5]
    lat = [67.5,70.25]
    time_step = 900
    duration = 1
    sep = 0.01
    start_time=[2022,6,10,12] #no data before 2022 5 21 00, seems that it changes from day to day, deleting old dates ?
    # start time particle advection [2022,6,10,12], [2022,6,12,12]
    #def __init__(self, lon=None, lat=None, time_step=None, sep=None, dur=None, date=None, hours=None, frames=5, jobs=2):
    '''
    ani = animations(time_step=0.5, sep=0.02, dur=15, frames=40)
    ani.LCS_and_particle_advection_gyre(3, 23, 'new_test')
    '''
    ani = animations(lon = lon, lat = lat, time_step = time_step, sep = sep, dur = duration, frames = 48)
    ani.LCS_and_particle_advection(start_time=start_time, outfile='small_test2')
    #LCS_and_particles_advection_animation('test', 48, lon, lat, start_time, time_step, duration, dxdy, LCS_start_time = 12)
