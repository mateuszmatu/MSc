import numpy as np
from cmath import inf
import matplotlib.pyplot as plt

def FTLE(file):
    d = np.load(file, allow_pickle='TRUE').item()
    x0 = d['initial_lon'][0,:]
    y0 = d['initial_lat'][:,0]
    duration = d['duration']
    dxdy = d['dxdy']*2

    X0, Y0 = np.meshgrid(x0,y0)
    
    LCS = {'lon':x0, 'lat':y0, 'RLCS': {}, 'ALCS':{}}

    for t in range(1, len(d['time'])+1):
        for fb, ar in zip(['f', 'b'], ['RLCS', 'ALCS']):
            x1, y1 = np.array(d['time'][t][fb]['lon']), np.array(d['time'][t][fb]['lat'])

            nx = X0.shape[0]
            ny = X0.shape[1]

            x1 = np.reshape(x1, (nx, ny))
            y1 = np.reshape(y1, (nx, ny))

            x = X0 - x1
            y = Y0 - y1

            dx = np.gradient(x)
            dy = np.gradient(y)

            F = np.zeros([nx,ny,2,2])
            F[:,:,0,0] = dx[0]/dxdy
            F[:,:,0,1] = dx[1]/dxdy
            F[:,:,1,0] = dy[0]/dxdy
            F[:,:,1,1] = dy[1]/dxdy

            largest_eig = np.zeros([nx,ny])

            for i in range(nx):
                for j in range(ny):
                    C = np.dot(np.transpose(F[i,j]), F[i,j])
                    eig = np.linalg.eigvals(C)
                    np.seterr(divide = 'ignore')
                    largest_eig[i,j] = np.log(np.max(eig))#np.log(np.sqrt(np.max(eig)))/duration
            largest_eig[largest_eig==-inf]=np.nan
            LCS[ar][t] = largest_eig
        LCS['ALCS'][t] = LCS['ALCS'][t][::-1,::-1]
    return LCS

if __name__ == '__main__':
    
    LCS = FTLE('flow_maps/20220319/values_m1.npy')
    LCS2 = FTLE('flow_maps/20220319/values_m2.npy')
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
    LCS = FTLE('flow_maps/double_gyre/values.npy')
    fig, ax = plt.subplots(2)
    a = ax[0].pcolormesh(LCS['lon'], LCS['lat'], LCS['RLCS'][1], cmap='jet')
    b = ax[1].pcolormesh(LCS['lon'], LCS['lat'], LCS['ALCS'][1], cmap='jet')

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
    