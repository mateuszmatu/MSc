from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs

def flow_map(file):
    transform = ccrs.NorthPolarStereo()
    plate = ccrs.PlateCarree()
    axs = plt.axes(projection = transform)

    for j in range(6):
        data = Dataset(f'{file}_m{j}.nc')
        lat, lon = data.variables['lat'], data.variables['lon']
        grid_size = len(lat)
        color = ['lightcoral', 'darkorange', 'olive', 'teal', 'violet', 
            'skyblue']
        plt.scatter(lon[:,1], lat[:,1], transform = plate, color=color[j], s = 25, label=f'Member {j}')
    plt.scatter(lon[:,0], lat[:,0], transform = plate, color='black', label='Initial position')
    plt.legend()
    axs.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
    axs.gridlines(draw_labels=True)
    axs.coastlines()
    plt.show()

if __name__ == '__main__':
    flow_map('flow_maps/20220323/f/20220323_f')