import cartopy.crs as ccrs
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy

def ftle(map, plot=False):
    data = Dataset(map)
    x0, y0, x1, y1, = data.variables['lat'][:,0], data.variables['lon'][:,0], data.variables['lat'][:,1], data.variables['lon'][:,1]
    size = len(x0)
    dx0 = 2*(x0[int(np.sqrt(size))]-x0[0])
    dy0 = 2*(y0[1]-y0[0])
    dx1, dy1 = np.array(np.gradient(x1)), np.array(np.gradient(y1))
    #something is wrong with dy1
    #for i in dy1:
    #    print(i)

    #Deformation-gradient tensor
    F = np.zeros([size, 2, 2])
    F[:, 0, 0] = dx1 / dx0
    F[:, 1, 1] = dy1 / dy0
    F[:, 0, 1] = dx1 / dy0
    F[:, 1, 0] = dy1 / dx0
    #Cauchy-Green strain tensor

    #remove the loop later on
    largest_eig = np.zeros(size)
    for i in range(size):
        C = np.dot(np.transpose(F[i]), F[i])
        eig = np.linalg.eigvals(C)
        largest_eig[i] = np.max(eig)
        
    largest_eig = largest_eig.reshape(int(np.sqrt(size)), int(np.sqrt(size)))

    if plot == True:
        transform = ccrs.NorthPolarStereo()
        plate = ccrs.PlateCarree()
        lons = np.linspace(y0[0]-0.1, y0[-1]+0.1, int(np.sqrt(size)))
        lats = np.linspace(x0[0]-0.1, x0[-1]+0.1, int(np.sqrt(size)))
        axs = plt.axes(projection = transform)
        plt.contourf(lons[1:-2], lats[1:-2], largest_eig[1:-2,1:-2], transform = plate, cmap = 'Reds')
        plt.colorbar()
        lat, lon = data.variables['lat'], data.variables['lon']
        axs.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
        axs.gridlines(draw_labels=True)
        axs.coastlines()
        plt.show()

if __name__ == '__main__':
    ftle('flow_maps/20220323/f/20220323_f_m0.nc', True)
