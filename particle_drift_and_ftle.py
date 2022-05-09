import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import cartopy
import cartopy.crs as ccrs

def particles_and_alcs(LCS, particles):
    d = xr.open_dataset(particles)
    lon = np.array(d['lon'])
    lat = np.array(d['lat'])
    transform = ccrs.NorthPolarStereo()
    plate = ccrs.PlateCarree()
    axs = plt.axes(projection = transform)
    plt.scatter(lon[:,-1], lat[:,-1], transform=plate)
    axs.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
    axs.gridlines(draw_labels=True)
    axs.coastlines()
    plt.show()


particles_and_alcs(1, 'flow_maps/20220419/b/particle_drift.nc')