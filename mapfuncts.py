from utils import *
from netCDF4 import Dataset, num2date

ensbounds = (-110.0, -46.0, -50.0, -11.0)
fig, ax = plt.subplots(1)
m = Basemap(ax = ax, projection='merc', llcrnrlon=ensbounds[0],
            urcrnrlon=ensbounds[1], llcrnrlat=ensbounds[2],
            urcrnrlat=ensbounds[3], resolution='l', area_thresh=5000)

fp = 'hgt.2015.nc'
data = Dataset('{:s}/Desktop/{:s}'.format(EV['HOME'],fp))
dates = num2date(data['time'][:], data['time'].units)
lat = data['lat'][:]
lon = data['lon'][:]
lons, lats = np.meshgrid(lon, lat)
lvls = data['level'][:]
heights = data['hgt'][:]

for lvl in [250, 500, 925, 1000]:
    hdex = lvls == lvl
    lines = m.contour(lons, lats, heights[0, hdex, :, :].squeeze(), latlon=True)
    clabels = plt.clabel(lines, fmt='%2.0f')
    fig.savefig(str(lvl))
