class initiate_particle_advection:
    def __init__(self, lon = None, lat = None, time_step = None, dxdy = 0.1, duration = None, date = None, hours = 2, start_time = 0):
        self.lon = lon
        self.lat = lat
        self.time_step = time_step
        self.dxdy = dxdy
        self.duration = duration
        self.date = date
        self.start_time = start_time

        x = np.arange(lon[0], lon[1], dxdy)
        y = np.arange(lat[0], lat[1], dxdy)
        self.lons, self.lats = np.meshgrid(x,y)

        self.path = f.create_directory(date)
        self.thredds_files = f.files_from_thredds(date)
        self.number_of_sims = int(hours/(self.duration/(60*60)))
    
    def set_up_dictionary(self, number_of_sims):
        #Creating a dictionary with values to be saved
        values = {'duration': self.duration, 'time_step': self.time_step, 'dxdy': self.dxdy, 'initial_lon': self.lons, 'initial_lat': self.lats, 'time': {}}
        for i in range(number_of_sims):
            values['time'][i+1] = {}
            values['time'][i+1]['f'] = {}
            values['time'][i+1]['b'] = {}
        return values

    @property
    def initiate_advection_all_members(self):
        #Run for all members for the current day
        for i in range(24):
            self.inititate_advection_one_member(i)

    def inititate_advection_one_member(self, member):
        #Run for a single member of the current day
        corr = f.correct_file(self.thredds_files, member)
        file = corr[0]
        _member = corr[1]
        
        values = self.set_up_dictionary(self.number_of_sims)
        for i in range(self.number_of_sims):
            for j in range(2):
                o = OceanDrift(loglevel=20)
                r = reader_netCDF_CF_generic.Reader(file, ensemble_member=_member)
                o.add_reader(r)
                if j == 0:
                    o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(hours=i*self.duration)+timedelta(hours=self.start_time))
                    o.run(duration=timedelta(seconds=self.duration), time_step=self.time_step, time_step_output=self.duration, outfile=f'{self.path}/m{member}_f.nc')
                if j == 1:
                    o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=r.start_time+timedelta(hours=i*self.duration)+timedelta(hours=self.start_time))
                    o.run(duration=timedelta(hours=self.duration), time_step=-self.time_step, time_step_output=-self.duration, outfile=f'{self.path}/m{member}_b.nc')
            _f = xr.open_dataset(f'{self.path}/m{member}_f.nc')
            _b = xr.open_dataset(f'{self.path}/m{member}_b.nc')
            os.remove(f'{self.path}/m{member}_f.nc')
            os.remove(f'{self.path}/m{member}_b.nc')
            values['time'][i+1]['f']['lon'] = np.array(_f['lon'][:,-1])
            values['time'][i+1]['f']['lat'] = np.array(_f['lat'][:,-1])
            values['time'][i+1]['b']['lon'] = np.array(_b['lon'][:,-1])
            values['time'][i+1]['b']['lat'] = np.array(_b['lat'][:,-1])
        np.save(f'{self.path}/values_m{member}.npy', values)
        
        return f'{self.path}/values_m{member}.npy'

    def aggregated(self, start_time, hours):
        file = 'https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_zdepth_be'
        start_time = datetime(start_time[0], start_time[1], start_time[2], start_time[3])

        duration = timedelta(hours=hours)
        values = {'duration': duration, 'time_step': self.time_step, 'dxdy': self.dxdy, 'initial_lon': self.lons, 'initial_lat': self.lats, 'time': {}}
        values['time'][1] = {}
        values['time'][1]['f'] = {}
        values['time'][1]['b'] = {}

        for j in range(2):
            o = OceanDrift(loglevel=20)
            r = reader_netCDF_CF_generic.Reader(file)
            o.add_reader(r)
            if j == 0:
                o.seed_elements(self.lons.ravel(), self.lats.ravel(), time = start_time)
                o.run(duration=duration, time_step=self.time_step, time_step_output=duration, outfile='f.nc')
            if j == 1:
                o.seed_elements(self.lons.ravel(), self.lats.ravel(), time = start_time)
                o.run(duration=duration, time_step=-self.time_step, time_step_output=duration, outfile='b.nc')
        _f = xr.open_dataset('f.nc')
        _b = xr.open_dataset('b.nc')
        os.remove('f.nc')
        os.remove('b.nc')
        values['time'][1]['f']['lon'] = np.array(_f['lon'][:,-1])
        values['time'][1]['f']['lat'] = np.array(_f['lat'][:,-1])
        values['time'][1]['b']['lon'] = np.array(_b['lon'][:,-1])
        values['time'][1]['b']['lat'] = np.array(_b['lat'][:,-1])
        np.save('values.npy', values)

    def aggregated_backwards(self, start_time, at_time = 0, duration=1, par=False, outfile='values'):
        file = 'https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_zdepth_be'
        start_time = datetime(start_time[0], start_time[1], start_time[2], start_time[3])

        o = OceanDrift(loglevel=20)
        r = reader_netCDF_CF_generic.Reader(file)
        o.add_reader(r)
        o.seed_elements(self.lons.ravel(), self.lats.ravel(), time=start_time+timedelta(hours = at_time))
        o.run(duration=timedelta(hours=duration), time_step=-self.time_step, time_step_output=timedelta(hours=duration), outfile=f'{outfile}.nc')

        _b = xr.open_dataset(f'{outfile}.nc')
        os.remove(f'{outfile}.nc')
        values = {'duration': duration, 'time_step': self.time_step, 'dxdy': self.dxdy, 'initial_lon': self.lons, 'initial_lat': self.lats, 
        'time': {1:{'f': {}, 'b':{}}}}
        values['time'][1]['b']['lon'] = np.array(_b['lon'][:,-1])
        values['time'][1]['b']['lat'] = np.array(_b['lat'][:,-1])
        values['time'][1]['f']['lon'] = np.zeros_like(_b['lon'][:,-1])
        values['time'][1]['f']['lat'] = np.zeros_like(_b['lat'][:,-1])

        np.save('values.npy', values)
        if par == True:
            path = f'tmp/{outfile}.npy'
            np.save(path, values)  
        else:
            path = f'{outfile}.npy'
            np.save(path)
        return path

        """  

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
                    largest_eig[i,j] = np.log(np.max(eig))#np.log(np.sqrt(np.max(eig)))/7200
            largest_eig[largest_eig==-inf]=np.nan
            LCS[ar][t] = largest_eig
        LCS['ALCS'][t] = LCS['ALCS'][t][::-1,::-1]
    return LCS
"""


def example_double_gyre(at_time = 3, duration = 15, time_step = 0.5, dxdy = 0.02, epsilon=0.25, omega=0.682, A=0.1):
    if not os.path.exists('flow_maps'):
        os.mkdir('flow_maps')
    if not os.path.exists('flow_maps/double_gyre'):
        os.mkdir('flow_maps/double_gyre')
    save_str = f'{at_time}_{duration}_{time_step}_{dxdy}_{epsilon}_{omega}_{A}'
    save_str = re.sub('\.', '', save_str)
    path_to_file = f'flow_maps/double_gyre/{save_str}.npy'    
    for i in range(2):
        o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('drift:advection_scheme', 'runge-kutta4')

        double_gyre = reader_double_gyre.Reader(epsilon=epsilon, omega=omega, A=A)
        print (double_gyre)

        o.add_reader(double_gyre)

        proj = double_gyre.proj
        x = np.arange(double_gyre.xmin, double_gyre.xmax, dxdy)
        y = np.arange(double_gyre.ymin, double_gyre.ymax, dxdy)

        X, Y = np.meshgrid(x,y)
        lons, lats = proj(X,Y, inverse=True)
        if i == 0:
            o.seed_elements(lons.ravel(), lats.ravel(),
                            time=double_gyre.initial_time+timedelta(seconds=at_time))
            o.run(duration=timedelta(seconds=duration), time_step=time_step, time_step_output=duration, outfile='flow_maps/double_gyre/double_gyre_f.nc')
        elif i == 1:
            o.seed_elements(lons.ravel(), lats.ravel(),
                            time=double_gyre.initial_time+timedelta(seconds=at_time))
            o.run(duration=timedelta(seconds=duration), time_step=-time_step, time_step_output=duration, outfile='flow_maps/double_gyre/double_gyre_b.nc')

    _f = xr.open_dataset('flow_maps/double_gyre/double_gyre_f.nc')
    _b = xr.open_dataset('flow_maps/double_gyre/double_gyre_b.nc')
    os.remove('flow_maps/double_gyre/double_gyre_f.nc')
    os.remove('flow_maps/double_gyre/double_gyre_b.nc')
    values = {'duration': duration, 'time_step': time_step, 'dxdy': dxdy, 'initial_lon': lons, 'initial_lat': lats, 
        'time': {1:{'f': {}, 'b':{}}}}
    values['time'][1]['f']['lon'] = np.array(_f['lon'][:,-1])
    values['time'][1]['f']['lat'] = np.array(_f['lat'][:,-1])
    values['time'][1]['b']['lon'] = np.array(_b['lon'][:,-1])
    values['time'][1]['b']['lat'] = np.array(_b['lat'][:,-1])
    np.save(f'flow_maps/double_gyre/{save_str}.npy', values)

    return path_to_file

def example_double_gyre_backward(at_time = 3, duration = 15, time_step = 0.5, dxdy = 0.02, epsilon=0.25, omega=0.682, A=0.1):
    if not os.path.exists('flow_maps'):
        os.mkdir('flow_maps')
    if not os.path.exists('flow_maps/double_gyre'):
        os.mkdir('flow_maps/double_gyre')
    save_str = f'{at_time}_{duration}_{time_step}_{dxdy}_{epsilon}_{omega}_{A}'
    save_str = re.sub('\.', '', save_str)
    path_to_file = f'flow_maps/double_gyre/b_{save_str}.npy' 

    o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
    o.set_config('environment:fallback:land_binary_mask', 0)
    o.set_config('drift:advection_scheme', 'runge-kutta4')

    double_gyre = reader_double_gyre.Reader(epsilon=epsilon, omega=omega, A=A)
    print (double_gyre)

    o.add_reader(double_gyre)

    proj = double_gyre.proj
    x = np.arange(double_gyre.xmin, double_gyre.xmax, dxdy)
    y = np.arange(double_gyre.ymin, double_gyre.ymax, dxdy)

    X, Y = np.meshgrid(x,y)
    lons, lats = proj(X,Y, inverse=True)
    o.seed_elements(lons.ravel(), lats.ravel(),
                    time=double_gyre.initial_time+timedelta(seconds=at_time))#timedelta(seconds=duration+at_time))
    o.run(duration=timedelta(seconds=duration), time_step=-time_step, time_step_output=duration, outfile='flow_maps/double_gyre/double_gyre_b.nc')

    _b = xr.open_dataset('flow_maps/double_gyre/double_gyre_b.nc')
    os.remove('flow_maps/double_gyre/double_gyre_b.nc')
    values = {'duration': duration, 'time_step': time_step, 'dxdy': dxdy, 'initial_lon': lons, 'initial_lat': lats, 
        'time': {1:{'f': {}, 'b':{}}}}
    values['time'][1]['b']['lon'] = np.array(_b['lon'][:,-1])
    values['time'][1]['b']['lat'] = np.array(_b['lat'][:,-1])
    values['time'][1]['f']['lon'] = np.zeros_like(_b['lon'][:,-1])
    values['time'][1]['f']['lat'] = np.zeros_like(_b['lat'][:,-1])

    np.save(f'flow_maps/double_gyre/b_{save_str}.npy', values)

    return path_to_file

def double_gyre_for_animation(at_time = 3, duration = 15, time_step = 0.5, dxdy = 0.02, epsilon=0.25, omega=0.682, A=0.1, outfile = None):
    if not os.path.exists('flow_maps'):
        os.mkdir('flow_maps')
    if not os.path.exists('flow_maps/double_gyre'):
        os.mkdir('flow_maps/double_gyre')
    path_to_file = f'flow_maps/double_gyre/{outfile}.npy' 

    o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
    o.set_config('environment:fallback:land_binary_mask', 0)
    o.set_config('drift:advection_scheme', 'runge-kutta4')

    double_gyre = reader_double_gyre.Reader(epsilon=epsilon, omega=omega, A=A)
    print (double_gyre)

    o.add_reader(double_gyre)

    proj = double_gyre.proj
    x = np.arange(double_gyre.xmin, double_gyre.xmax, dxdy)
    y = np.arange(double_gyre.ymin, double_gyre.ymax, dxdy)

    X, Y = np.meshgrid(x,y)
    lons, lats = proj(X,Y, inverse=True)
    o.seed_elements(lons.ravel(), lats.ravel(),
                    time=double_gyre.initial_time+timedelta(seconds=at_time))#timedelta(seconds=duration+at_time))
    o.run(duration=timedelta(seconds=duration), time_step=-time_step, time_step_output=duration, outfile=f'flow_maps/double_gyre/{outfile}.nc')

    _b = xr.open_dataset(f'flow_maps/double_gyre/{outfile}.nc')
    os.remove(f'flow_maps/double_gyre/{outfile}.nc')
    values = {'duration': duration, 'time_step': time_step, 'dxdy': dxdy, 'initial_lon': lons, 'initial_lat': lats, 
        'time': {1:{'f': {}, 'b':{}}}}
    values['time'][1]['b']['lon'] = np.array(_b['lon'][:,-1])
    values['time'][1]['b']['lat'] = np.array(_b['lat'][:,-1])
    values['time'][1]['f']['lon'] = np.zeros_like(_b['lon'][:,-1])
    values['time'][1]['f']['lat'] = np.zeros_like(_b['lat'][:,-1])

    np.save(f'flow_maps/double_gyre/{outfile}.npy', values)

    return path_to_file

    def LCS_and_particles_advection_animation(outfile, frames, lon, lat, start_time, time_step, duration, dxdy, LCS_start_time = 0):
    '''
        Creates animation of an area of advected particles and the LCS field.
    Args:
        outfile [str]           : a string of the outfile name.
        frames [int]            : number of frames in animation.
        lon [list]              : list of min lon and max lon.
        lat [list]              : list of min lat and max lat.
        start_time [list]       : list of start time [year, month, day, hour]
        time_step [int]         : the time step in seconds of LCS simulation
        duration [int]          : the duration of LCS simulation in hours
        dxdy [float]            : separation between particles in LCS
        LCS_start_time [int]    : start time of LCS field compared to particle simulation in hours.
    '''

    def par(i):
        ip = ipa(lon, lat, time_step, dxdy, duration)
        file = ip.aggregated_backwards(start_time, at_time=i, par=True, outfile=i)
        LCS = FTLE(file)
        return LCS['ALCS'][1]
    
    def animation(i):
        x = []
        x = arr[i]
        posx = []
        posy = []
        posx = pos['lon'][:,i+LCS_start_time]
        posy = pos['lat'][:,i+LCS_start_time]
        ax.set_title(f'frame {i}/{frames}')
        ax.pcolormesh(in_x, in_y, x, cmap='jet', vmin=-0.2, vmax=0.5, transform=plate)
        ax.scatter(posx, posy, s=3, color='black', transform=plate)
        ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
        ax.coastlines()

    transform = ccrs.NorthPolarStereo()
    plate = ccrs.PlateCarree()
    fig = plt.figure(figsize=(20,20))
    ax = plt.axes(projection=transform)
    pos = xr.load_dataset('particle_drift.nc')

    arr = Parallel(n_jobs=2)(delayed(par)(i) for i in range(frames))
    data = np.load('tmp/1.npy',allow_pickle='TRUE').item()
    in_x = data['initial_lon']
    in_y = data['initial_lat']
    ani=FuncAnimation(fig, func=animation, interval=300,frames=frames)
    ani.save(f'gifs/{outfile}.gif', dpi=80,writer='imagemagick')

    def remove_files(i):
        os.remove(f'tmp/{i}.npy')

    Parallel(n_jobs=4)(delayed(remove_files)(i) for i in range(frames))

def LCS_and_particles_advection_animation_gyre(outfile, frames, t0, t1, duration = 15, time_step = 0.5, dxdy = 0.02, epsilon=0.25, omega=0.682, A=0.1):

    def animation(i):
        LCS = arr[i]
        x = []
        x = LCS
        posx = []
        posy = []
        posx = pos['lon'][:,i]
        posy = pos['lat'][:,i]
        ax.clear()
        ax.set_title(f'frame {i}/{frames}')
        ax.pcolormesh(in_x, in_y, x, cmap='jet')
        ax.scatter(posx, posy, s=5, color='black')

    def advect_particles():
        o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
        o.set_config('environment:fallback:land_binary_mask', 0)
        o.set_config('drift:advection_scheme', 'runge-kutta4')

        double_gyre = reader_double_gyre.Reader(epsilon=epsilon, omega=omega, A=A)
        print (double_gyre)

        o.add_reader(double_gyre)

        x = [.9]
        y = [.5]
        lon, lat = double_gyre.xy2lonlat(x, y)
        time_step = times[1]-times[0]
        o.seed_elements(lon, lat, radius=.1, number=5000,
                        time=double_gyre.initial_time+timedelta(seconds=t0))
        o.run(duration=timedelta(seconds=t1), time_step=time_step, outfile='flow_maps/double_gyre/particle_advection.nc')

    def par(i):
        file = dgfa(times[i],duration,time_step,dxdy,epsilon,omega,A,outfile=i)
        LCS = FTLE(file)
        return LCS['ALCS'][1]
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    times = np.linspace(t0,t1,frames)
    advect_particles()
    pos = xr.load_dataset('flow_maps/double_gyre/particle_advection.nc')

    arr = Parallel(n_jobs=4)(delayed(par)(i) for i in range(frames))
    data = np.load('flow_maps/double_gyre/1.npy',allow_pickle='TRUE').item()
    in_x = data['initial_lon']
    in_y = data['initial_lat']
    ani = FuncAnimation(fig, func=animation, interval=200, frames=frames)
    ani.save(f'gifs/{outfile}.gif', dpi = 80, writer = 'imagemagick')

    def remove_files(i):
        os.remove(f'flow_maps/double_gyre/{i}.npy')

    Parallel(n_jobs=4)(delayed(remove_files)(i) for i in range(frames))