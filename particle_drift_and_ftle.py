from example_double_gyre import double_gyre_for_animation as dgfa
from grid_advection import initiate_particle_advection as ipa
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

def particles_and_alcs(LCS, particles):
    d = xr.open_dataset(particles)
    lon = np.array(d['lon'])
    lat = np.array(d['lat'])
    transform = ccrs.NorthPolarStereo()
    plate = ccrs.PlateCarree()
    plt.figure(figsize=(20,20))
    axs = plt.axes(projection = transform)
    plt.pcolormesh(LCS['lon'], LCS['lat'], LCS['ALCS'][1], cmap='jet', transform=plate, vmin=-0.2, vmax=0.5)
    plt.colorbar(shrink=0.7)

    plt.scatter(lon[:,-1], lat[:,-1], transform=plate, color='black', s=3)
    axs.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
    axs.gridlines(draw_labels=True)
    axs.coastlines()
    plt.savefig('plot.png')

def lcs_plot(LCS, outfile='lcsplot'):
    transform = ccrs.NorthPolarStereo()
    plate = ccrs.PlateCarree()
    plt.figure(figsize=(20,20))
    axs = plt.axes(projection=transform)
    axs.pcolormesh(LCS['lon'], LCS['lat'], LCS['ALCS'][1], cmap='jet', transform=plate, vmin=-0.2, vmax=0.5)

    axs.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
    axs.gridlines(draw_labels=True)
    axs.coastlines()
    plt.savefig(f'{outfile}.png')

def LCS_and_particles_advection_animation(outfile, _frames, _start_time):

    def par(i):
        lon = [12.5,16]
        lat = [68,70]
        time_step = 900
        duration = 3600
        dxdy = 0.02
        start_time = 60*60*24*2
        ip = ipa(lon,lat,time_step,dxdy,duration,'20220419',hours=2,start_time=start_time)
        file=ip.aggregated_backwards(_start_time,at_time=i, par=True,outfile=i)
        LCS = FTLE(file)
        return LCS['ALCS'][1]
    
    def animation(i):
        x = []
        x = arr[i]
        posx = []
        posy = []
        posx = pos['lon'][:,i]
        posy = pos['lat'][:,i]
        ax.set_title(f'frame {i}/{_frames}')
        ax.pcolormesh(in_x, in_y, x, cmap='jet', vmin=-0.2, vmax=0.5, transform=plate)
        ax.scatter(posx, posy, s=3, color='black', transform=plate)
        ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
        ax.coastlines()

    transform = ccrs.NorthPolarStereo()
    plate = ccrs.PlateCarree()
    fig = plt.figure(figsize=(20,20))
    ax = plt.axes(projection=transform)
    pos = xr.load_dataset('particle_drift.nc')

    arr = Parallel(n_jobs=2)(delayed(par)(i) for i in range(_frames))
    data = np.load('tmp/1.npy',allow_pickle='TRUE').item()
    in_x = data['initial_lon']
    in_y = data['initial_lat']
    ani=FuncAnimation(fig, func=animation, interval=200,frames=_frames)
    ani.save(f'gifs/{outfile}.gif', dpi=80,writer='imagemagick')

    

    """
    def animation(i):
        lon = [7.5,17]
        lat = [68,72.3]
        time_step = 900
        duration = 7200
        dxdy = 0.01
        start_time = 60*60*24*2
        ip = ipa(lon,lat,time_step,dxdy,duration,'20220419', hours=2, start_time=start_time)
        file = ip.aggregated_backwards(_start_time,i)
        LCS = FTLE(file)
        x = []
        x = LCS['ALCS'][1]
        posx = []
        posy = []
        posx = pos['lon'][:,i]
        posy = pos['lat'][:,i]
        ax.clear()
        #ax.set_title(f'frame {i}/{frames}')
        ax.pcolormesh(LCS['lon'], LCS['lat'], x, cmap='jet', vmin=-0.2, vmax=0.5, transform=plate)
        ax.scatter(posx, posy, s=3, color='black', transform=plate)
        ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
        ax.coastlines()

    transform = ccrs.NorthPolarStereo()
    plate = ccrs.PlateCarree()
    fig = plt.figure(figsize=(20,20))
    ax = plt.axes(projection=transform)
    pos = xr.load_dataset('particle_and_lcs/particle_drift_two_weeks.nc')

    ani = FuncAnimation(fig, func=animation, interval=200, frames=_frames)
    ani.save(f'gifs/{outfile}.gif', dpi = 80, writer='imagemagick')
    """

def LCS_eps(lon, lat, time_step, dxdy, duration, date, hours, start_time, members):
    for member in members:
        ip = ipa(lon,lat,time_step,dxdy,duration,date,hours,start_time)
        file = ip.inititate_advection_one_member(member)
        LCS = FTLE(file)
        lcs_plot(LCS, f'flow_maps/{date}/lcs_plot_m{member}')


"""
LCS = FTLE('values.npy')
particles_and_alcs(LCS, 'particle_and_lcs/particle_drift_two_weeks.nc')

"""
LCS_and_particles_advection_animation(outfile='test', _frames=10, _start_time=[2022,5,20,12])

#LCS_eps([7.5,16.5], [68,72],900, 0.01, 7200,'20220530', 2, 4, [0,1,2,3,4,5])
"""

LCS = FTLE('flow_maps/20220419/values_m1.npy')
LCS2 = FTLE('flow_maps/20220419/values_m2.npy')
fig, ax = plt.subplots(2,2)
a = ax[0,0].pcolormesh(LCS['lon'], LCS['lat'], LCS['RLCS'][1], cmap='jet')
c = ax[1,0].pcolormesh(LCS['lon'], LCS['lat'], LCS['ALCS'][1], cmap='jet', vmin=-1, vmax=2)
b = ax[0,1].pcolormesh(LCS2['lon'], LCS['lat'], LCS['RLCS'][1], cmap='jet')
d = ax[1,1].pcolormesh(LCS2['lon'], LCS['lat'], LCS['ALCS'][1], cmap='jet',vmin=-1, vmax=2)

ax[0,0].set(title='Repelling LCS')
ax[0,1].set(title='Repelling LCS')
ax[1,0].set(title='Attracting LCS')
ax[1,1].set(title='Attracting LCS')
l = np.array([[a,b],[c,d]])
from mpl_toolkits.axes_grid1 import make_axes_locatable
for j in range(2):
    for i in range(2):
        divider = make_axes_locatable(ax[i,j])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(l[i,j], cax=cax)
plt.show()
"""