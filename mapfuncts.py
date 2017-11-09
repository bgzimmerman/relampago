from utils import *
from netCDF4 import Dataset, num2date
from scipy.ndimage.filters import gaussian_filter

ensbounds = (-110.0, -46.0, -50.0, -11.0) #lon,lon,lat,lat
def gen_basemap(ax):
    m = Basemap(ax = ax, projection='merc', llcrnrlon=ensbounds[0],
                urcrnrlon=ensbounds[1], llcrnrlat=ensbounds[2],
                urcrnrlat=ensbounds[3], resolution='l', area_thresh=5000)
    m.drawcoastlines(linewidth=2)
    m.drawcountries(linewidth=2)
    m.drawstates(linewidth=1, color='0.8')
    parallels = np.arange(-60,-5, 5)
    meridians = np.arange(-120,-30, 10)
    m.drawparallels(parallels,labels=[False,True,True,False])
    m.drawmeridians(meridians,labels=[True,False,True,False])
    return ax, m

fp = 'hgt.2015.nc'
data = Dataset('{:s}/Desktop/{:s}'.format(EV['HOME'],fp))
dates = num2date(data['time'][:], data['time'].units)

lat = data['lat'][:]
lon = data['lon'][:]
lons, lats = np.meshgrid(lon, lat)
lvls = data['level'][:]
heights = data['hgt'][:]

pltdat = {
        '250' : np.arange(8000, 13000, 120),
        '500' : np.arange(4000, 7000, 60),
        '925' : np.arange(0, 4000, 30),
        '1000': np.arange(0, 4000, 20)
        }


brief_date = datetime(2015, 11, 12, 11, 0, 0)
tdex = 0

times = [date for date in dates if (date.month == 11 and (date.day == 11 or date.day == 12))]

# wind data
fp = 'uwnd.2015.nc'
data = Dataset('{:s}/Desktop/{:s}'.format(EV['HOME'],fp))
uwnd = data['uwnd']

fp = 'vwnd.2015.nc'
data = Dataset('{:s}/Desktop/{:s}'.format(EV['HOME'],fp))
vwnd = data['vwnd']

def roundup(x):
    return int(np.ceil(x/10.0)) * 10

for instance in times:
    tdex = dates == instance
    for lvl in [250, 500, 925, 1000]:
        fig, ax = plt.subplots(1, figsize = (12,10))
        ax, m = gen_basemap(ax)
        hdex = lvls == lvl
        data = gaussian_filter(heights[tdex, hdex, :, :].squeeze(), 1)
        lines = m.contour(lons, lats, data,
                            pltdat[str(lvl)], colors = 'k', latlon=True)
        clabels = plt.clabel(lines, fmt='%2.0f')

        udata = gaussian_filter(uwnd[tdex, hdex, :, :].squeeze(), 1)
        vdata = gaussian_filter(vwnd[tdex, hdex, :, :].squeeze(), 1)

        speed = np.sqrt(udata**2 + vdata**2)

        cs = m.contourf(lons, lats, speed, levels = np.arange(0, roundup(speed.max()), 5), latlon=True)
        m.barbs(lons, lats, udata, vdata, latlon = True)
        cb = plt.colorbar(cs, orientation = 'horizontal')
        cb.set_label('$m\/s^{-1}$', size = 'x-large')
        fig.suptitle('{:%Y%m%d%H}_{:d}hPa Heights and Winds'.format(instance, lvl))
        fig.savefig('{:%Y%m%d%H}_{:d}'.format(instance, lvl))
        plt.close(fig)


