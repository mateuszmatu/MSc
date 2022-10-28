from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.oceandrift import OceanDrift
import numpy as np
from datetime import timedelta
from functions import create_directory, files_from_thredds, correct_file
from datetime import datetime
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def initiate_files_all_members(lon, lat, time_step, dxdy, duration, date):
    for i in range(24):
        initiate_files_one_member(lon,lat,time_step,dxdy,duration,date,i)

def initiate_files_one_member(lon, lat, time_step, dxdy, duration, date, member):
    path = create_directory(date)
    thredds_files = files_from_thredds(date)
    corr = correct_file(thredds_files, member)
    file = corr[0]
    _member = corr[1]
    
    o = OceanDrift(loglevel=20)
    r = reader_netCDF_CF_generic.Reader(file, ensemble_member = _member)
    o.add_reader(r)
    x = np.arange(lon[0], lon[1], dxdy)
    y = np.arange(lat[0], lat[1], dxdy)
    #proj = r.proj
    lons, lats = np.meshgrid(x,y)
    #lons, lats = proj(X,Y, inverse=True)
    o.seed_elements(lons.ravel(), lats.ravel(), time=r.start_time)
    o.run(duration=timedelta(seconds=duration), time_step=time_step, outfile=f'{path}/particle_drift_m{member}.nc')

def advection_from_aggregated_file(lon, lat, time_step, dxdy, start_time = None, end_time = None):
    '''
        Advects particles using the aggregated barents eps zdepth file
    Args:
        lon [list] : 
        lat [list] : 
        time_step [int] : 
        dxdy 
    '''
    file = 'https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_zdepth_be'

    o = OceanDrift(loglevel=30)
    r = reader_netCDF_CF_generic.Reader(file)
    o.add_reader(r)
    x = np.arange(lon[0], lon[1], dxdy)
    y = np.arange(lat[0], lat[1], dxdy/2)
    lons, lats = np.meshgrid(x,y)

    if start_time is None:
        start_time = r.start_time
    else:
        start_time = datetime(start_time[0], start_time[1], start_time[2], start_time[3])
    
    if end_time is None:
        end_time = start_time + timedelta(hours = 72)
    else:
        end_time = datetime(end_time[0], end_time[1], end_time[2], end_time[3])
    duration = end_time - start_time
    o.seed_elements(lons.ravel(), lats.ravel(), time = start_time)
    o.run(duration = duration, time_step=time_step, outfile='particle_drift.nc')


def velocity(lon, lat, date, member, step=2):

    path = create_directory(date)
    thredds_files = files_from_thredds(date)
    corr = correct_file(thredds_files, member)
    file = corr[0]
    _member = corr[1]

    data = xr.open_dataset(file).isel(ensemble_member=_member)

    a1 =100
    a2 =300
    a3 = 400
    a4 = 550
    u = np.array(data['u'][12,0,a1:a2:step,a3:a4:step])
    v = np.array(data['v'][12,0,a1:a2:step,a3:a4:step])
    lons = np.array(data['lon'][a1:a2:step,a3:a4:step])
    lats = np.array(data['lat'][a1:a2:step,a3:a4:step])
    
    #lon = np.linspace(lon[0], lon[1], len(u))
    #lat = np.linspace(lat[0], lat[1], len(v))
    #lons, lats = np.meshgrid(lon, lat)
    plt.figure(figsize=(20,20))
    ax = plt.axes(projection = ccrs.NorthPolarStereo())
    plt.quiver(lons, lats, u, v, transform = ccrs.PlateCarree(), scale = 25)
    ax.coastlines()
    ax.gridlines(draw_labels=True)
    plt.title(f'Velocities at Vesterålen member {member} \n 30.05.22 12.00', fontdict={'fontsize': 30})
    plt.savefig(f'{path}/vel_m{member}.png', bbox_inches='tight')


    

if __name__ == '__main__':
    lon = [15,17]
    lat = [69.25,70]

    time_step = 3600
    dxdy = 0.01
    #initiate_file_all_members(lon, lat, time_step, dxdy, duration, '20220420')
    #initiate_files_one_member(lon, lat, time_step, dxdy, duration, '20220419', 0)
    advection_from_aggregated_file(lon, lat, time_step, dxdy, [2022,8,1,12], [2022,8,7,12])
    """
    for i in range(6):
        velocity(lon,lat,'20220530', i)
    """