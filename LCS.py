from example_double_gyre import example_double_gyre
import numpy as np
from cmath import inf
import matplotlib.pyplot as plt
from grid_advection import particle_grid_displacement as pgd

def FTLE(file):
    '''
        A script for calculating attracting and repelling LCS fields using the Couchy-Green Strain-Tensor and FTLE method.
    Args:
        file    [str]   :   Filename
    Returns:
        LCS     [dict]  :   A dictionary containing initial lons and lats, and calculated RLCS and ALCS
    '''

    d = np.load(file, allow_pickle='TRUE').item()
    x0 = d['in_lons'][0,:]
    y0 = d['in_lats'][:,0]
    duration = d['duration']
    sep = d['in_separation']*2

    X0, Y0 = np.meshgrid(x0,y0)
    
    LCS = {'lon':x0, 'lat':y0}

    for fb, ar in zip(['f', 'b'], ['RLCS', 'ALCS']):
        x1, y1 = np.array(d[fb]['lon']), np.array(d[fb]['lat'])
        nx, ny = X0.shape[0], X0.shape[1]
        x1, y1 = np.reshape(x1, (nx, ny)), np.reshape(y1, (nx, ny))
        x, y = X0 - x1, Y0 - y1

        dx, dy = np.gradient(x), np.gradient(y)

        F = np.zeros([nx,ny,2,2])
        F[:,:,0,0] = dx[0]/sep
        F[:,:,0,1] = dx[1]/sep
        F[:,:,1,0] = dy[0]/sep
        F[:,:,1,1] = dy[1]/sep

        largest_eig = np.zeros([nx,ny])

        for i in range(nx):
            for j in range(ny):
                C = np.dot(np.transpose(F[i,j]), F[i,j])
                eig = np.linalg.eigvals(C)
                np.seterr(divide = 'ignore')
                largest_eig[i,j] = np.log(np.max(eig))
        largest_eig[largest_eig==-inf]=np.nan
        LCS[ar] = largest_eig
    LCS['ALCS'] = LCS['ALCS'][::-1,::-1]

    return LCS
    
if __name__ == '__main__':
    
    import cartopy
    import cartopy.crs as ccrs

    r = pgd([14,17.5], [69,70.25], 900, 0.01, 1, '20220610', 12)
    files = r.displacement_all_members(forwards = False, num_of_members=24)
    
    transform = ccrs.NorthPolarStereo()
    plate = ccrs.PlateCarree()
    fig, ax = plt.subplots(3,2,figsize=(10,10),subplot_kw={'projection': transform})
    LCS = []
    for file in files:
        LCS.append(FTLE(file))

    member = 0
    for i in range(2):
        for j in range(2):
            ax[i,j].pcolormesh(LCS[member]['lon'], LCS[member]['lat'], LCS[member]['ALCS'], cmap='jet', vmin=-0.2, vmax=0.35, transform=plate)
            ax[i,j].add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
            ax[i,j].coastlines()
            ax[i,j].gridlines(draw_labels=True)
            member+=1

    sum24 = 0
    sum4 = 0
    for i in range(len(LCS)):
        sum24 += LCS[i]['ALCS']
    mean24 = sum24/len(LCS)
    for i in range(4):
        sum4 += LCS[i]['ALCS']
    mean4 = sum4/4
    ax[2,0].pcolormesh(LCS[0]['lon'], LCS[0]['lat'], mean4, cmap='jet', vmin=-0.2, vmax=0.35, transform=plate)
    ax[2,0].add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
    ax[2,0].coastlines()
    ax[2,0].gridlines(draw_labels=True)
    ax[2,1].pcolormesh(LCS[0]['lon'], LCS[0]['lat'], mean24, cmap='jet', vmin=-0.2, vmax=0.35, transform=plate)
    ax[2,1].add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
    ax[2,1].coastlines()
    ax[2,1].gridlines(draw_labels=True)
    plt.tight_layout()
    plt.show()

    """
    file = example_double_gyre(forwards=False, sep=0.005)
    LCS = FTLE(file)
    plt.pcolormesh(LCS['lon'], LCS['lat'], LCS['ALCS'], cmap='plasma')
    plt.colorbar()
    plt.show()
    """
    """
    LCS = FTLE('flow_maps/20220529/lon12-5_16_lat68_70_h12_m0.npy')
    LCS2 = FTLE('flow_maps/20220529/lon12-5_16_lat68_70_h12_m1.npy')
    fig, ax = plt.subplots(2,2)
    a = ax[0,0].pcolormesh(LCS['lon'], LCS['lat'], LCS['RLCS'], cmap='ocean')
    c = ax[1,0].pcolormesh(LCS['lon'], LCS['lat'], LCS['ALCS'], cmap='gist_earth', vmin=0, vmax=2)
    b = ax[0,1].pcolormesh(LCS2['lon'], LCS2['lat'], LCS2['RLCS'], cmap='ocean')
    d = ax[1,1].pcolormesh(LCS2['lon'], LCS2['lat'], LCS2['ALCS'], cmap='gist_earth',vmin=0, vmax=2)

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
    