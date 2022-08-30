from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import xarray as xr
import os
import numpy as np
import matplotlib.cm as cm

def flow_map_all_members(path, fb='f'):
    transform = ccrs.NorthPolarStereo()
    plate = ccrs.PlateCarree()
    axs = plt.axes(projection = transform)
    files = []
    for file in os.listdir(path):
        files.append(file)
    files.sort()
    colors = cm.rainbow(np.linspace(0, 1, 6))
    m=0
    for file, c in zip(files, colors):
        data = np.load(f'{path}/{file}', allow_pickle='TRUE').item()
        times = len(data['time'])
        for time in range(1, 2):
            lon = data['time'][time][fb]['lon']
            lat = data['time'][time][fb]['lat']
            plt.scatter(lon, lat, transform=plate, s=25, color=c, label=m)
        m+=1
    in_lon = data['initial_lon']
    in_lat = data['initial_lat']
    plt.scatter(in_lon, in_lat, transform=plate, color='black', label='Initial positions')
    axs.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
    plt.legend()
    axs.gridlines(draw_labels=True)
    axs.coastlines()
    plt.show()

def flow_map_one_member(file):
    transform = ccrs.NorthPolarStereo()
    plate = ccrs.PlateCarree()
    axs = plt.axes(projection = transform)

    data = Dataset(file)
    lat, lon = data.variables['lat'], data.variables['lon']
    grid_size = len(lat)
    plt.scatter(lon[:,0], lat[:,0], transform = plate, color='black', label='Initial position')
    plt.scatter(lon[:,1], lat[:,1], transform = plate, s = 25)
    
    plt.legend()
    axs.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
    axs.gridlines(draw_labels=True)
    axs.coastlines()
    plt.show()

def flow_map_double_gyre(file):
    data = xr.open_dataset(file)
    lat, lon = data.variables['lat'], data.variables['lon']

    return lon[:,0], lat[:,0], lon[:,1], lat[:,1]

def flow_map_selected_members(path, members = [0]):
    transform = ccrs.NorthPolarStereo()
    plate = ccrs.PlateCarree()
    axs = plt.axes(projection = transform)
    files = []
    for member in members:
        files.append(f'{path}values_m{member}.npy')
    
    files.sort()
    colors = cm.rainbow(np.linspace(0, 1, len(members)))
    m=0
    for file, c in zip(files, colors):
        data = np.load(f'{file}', allow_pickle='TRUE').item()
        for time in range(1, 2):
            lon = data['time'][1]['b']['lon']
            lat = data['time'][1]['b']['lat']
            plt.scatter(lon, lat, transform=plate, s=25, color=c, label=m)
        m+=1
    in_lon = data['initial_lon']
    in_lat = data['initial_lat']
    plt.scatter(in_lon, in_lat, transform=plate, color='black', label='Initial positions')
    axs.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
    plt.legend(bbox_to_anchor=(1,1), loc="upper left")
    axs.gridlines(draw_labels=True)
    axs.coastlines()
    plt.savefig('flow_map.png')
    
if __name__ == '__main__':
    #flow_map_all_members('flow_maps/20220323/f/20220323_f')
    flow_map_selected_members('flow_maps/20220529/', [0,1,2,3,4,5])