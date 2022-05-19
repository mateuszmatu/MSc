
from cmath import inf
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle
import xarray as xr



def FTLE(map, initial_values):

    data = xr.open_dataset(map)
    initial_values = np.load(initial_values,allow_pickle='TRUE').item()
    x0 = initial_values['initial_lons'][0,:]
    y0 = initial_values['initial_lats'][:,0]
    duration = initial_values['duration']
    dxdy = initial_values['dxdy']*2

    X0, Y0 = np.meshgrid(x0,y0)
    x1, y1 = np.array(data['lon'][:,1]), np.array(data['lat'][:,1])

    nx = X0.shape[0]
    ny = X0.shape[1]
    x1 = np.reshape(x1, (nx, ny))
    y1 = np.reshape(y1, (nx, ny))

    x = X0-x1
    y = Y0-y1
    
    dx = np.gradient(x)
    dy = np.gradient(y)

    F = np.zeros([nx,ny,2,2])
    F[:,:,0,0] = dx[0]/dxdy
    F[:,:,0,1] = dx[1]/dxdy
    F[:,:,1,0] = dy[0]/dxdy
    F[:,:,1,1] = dy[1]/dxdy

    largest_eig = np.zeros([nx,ny])
   
    for i in tqdm(range(nx)):
        for j in range(ny):
            C = np.dot(np.transpose(F[i,j]), F[i,j])
            eig = np.linalg.eigvals(C)
            np.seterr(divide = 'ignore')
            largest_eig[i,j] = np.log(np.max(eig))#np.log(np.sqrt(np.max(eig)))/duration

    largest_eig[largest_eig==-inf]=np.nan
    ftle={'lcs':largest_eig, 'lons':x0, 'lats':y0}
    return ftle

if __name__ == '__main__':
    
    vals = 'flow_maps/double_gyre/initial_values.npy'
    rep = FTLE('flow_maps/double_gyre/double_gyre_f.nc', vals)
    att = FTLE('flow_maps/double_gyre/double_gyre_b.nc', vals)
    #attracting = attracting[::-1,::-1] #flipping for some reason

    fig, ax = plt.subplots(2)
    a = ax[0].pcolormesh(rep['lons'], rep['lats'], rep['lcs'], cmap='jet')
    b = ax[1].pcolormesh(att['lons'], att['lats'], att['lcs'], cmap='jet')
    ax[0].set(title='Repelling LCS')
    ax[1].set(title='Attracting LCS')
    l = [a,b]
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    for i in range(2):
        divider = make_axes_locatable(ax[i])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(l[i], cax=cax)
    plt.show()
    """
    vals = 'flow_maps/20220419/initial_values.npy'
    repelling = FTLE('flow_maps/20220419/f/m1.nc', vals)
    attracting = FTLE('flow_maps/20220419/b/m1.nc', vals)
    #attracting = attracting[::-1,::-1] #flipping for some reason

    fig, ax = plt.subplots(2)
    a = ax[0].imshow(repelling['lcs'], interpolation='nearest', origin='lower', cmap='jet')
    b = ax[1].imshow(attracting['lcs'], interpolation='nearest', origin='lower', cmap='jet')
    ax[0].set(title='Repelling LCS')
    ax[1].set(title='Attracting LCS')
    l = [a,b]
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    for i in range(2):
        divider = make_axes_locatable(ax[i])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(l[i], cax=cax)
    plt.show()
    """