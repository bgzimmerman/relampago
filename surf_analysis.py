import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import cm
import os
from mpl_toolkits.basemap import Basemap
EV = os.environ
# The index in the downloaded file is in spanish
translate = True
english_columns = ['lat', 'lon', 'station', 'UTC', 'temp', 'humidity',
                'dewpoint', 'windspeed', 'winddir', 'prcp24', 'prcp6',
                'prcp1', 'pressure']

#Set up path to data files
dir_path = '{:s}/relampago/AR_surf_data'.format(EV['HOME'])
testfile = dir_path + '/2017101916.csv'

# Load the data
data = pd.read_csv(testfile)

# Translate the columns
data.columns = english_columns if translate else data.columns

# Get the set of unique stations - some have data every 5 minutes
stations = set(data['station'])

# Group the data by stations so we can perform averaging
rawdata = data.groupby('station')

# The data file is for every hour, so we want to average the # observations within the time period, and fill a new dataframe
hourlydata = pd.DataFrame(columns = stations)

# Set the variables we want to evaluate

variables = ['lat', 'lon', 'temp', 'humidity', 'dewpoint',
                'windspeed', 'winddir', 'prcp24', 'prcp6',
                'prcp1', 'pressure']

for station in stations:
    print station
    for var in variables:
        hourlydata[station][var] = rawdata[var].mean()[station]

# Now need to create some grids and fill them in with the data

bnds = (-43.0, -20.0, -73.0, -50.0) # (minlat, maxlat, minlon, maxlon)
step = 0.2
grid_lat, grid_lon = np.mgrid[bnds[0]:bnds[1]:step, bnds[2]:bnds[3]:step]

grids = {}
# This is the interpolation part of the code
for var in variables:
    points = np.zeros((len(stations), 2))
    values = np.zeros(len(stations))
    for i, station in enumerate(stations):
        points[i, 0] = hourlydata[station]['lat']
        points[i, 1] = hourlydata[station]['lon']
        values[i] = hourlydata[station][var]

        # Get rid of nan values
    nan_index = np.isnan(values) | (values <= 0)
    if nan_index.sum() == 398:
        continue
    points = points[~nan_index,:]
    values = values[~nan_index]

    grids[var] = griddata(points, values, (grid_lat, grid_lon), method = 'linear')

    ensbounds=(-72.0, -48.0, -42.0, -18.0)
    fig, ax = plt.subplots(1)
    m = Basemap(projection='merc', llcrnrlon=ensbounds[0],
                urcrnrlon=ensbounds[1], llcrnrlat=ensbounds[2],
                urcrnrlat=ensbounds[3], resolution='l', area_thresh=5000)
    m.drawcoastlines()
    m.drawcountries()
    parallels = np.arange(-40., -19., 5.)
    meridians = np.arange(-70., -49., 5.)
    m.drawparallels(parallels, labels = [True, False, False, False])
    m.drawmeridians(meridians, labels = [False, False, False, True])


    cs = m.pcolormesh(grid_lon, grid_lat, grids[var], latlon = True, cmap = cm.RdBu_r)
    x, y = m(points[:,1], points[:,0])


    cormenlat = (-31.42,-32.89)
    cormenlon = (-64.18,-68.84)
    m.scatter(x, y, s = 1, marker = 'x', color = 'k', alpha = 0.8)
    xpt, ypt = m(cormenlon, cormenlat)
    m.scatter(xpt, ypt, s = 20, marker = 'o', color = 'r')
    def get_terrain(terrfile='/home/disk/meso-home/bzim/relampago/geo_em.d01.nc'):
        """
        Loads terrain height from WRF file
        """
        from scipy.ndimage.filters import gaussian_filter
        from netCDF4 import Dataset
        with Dataset(terrfile,'r') as dset:
            lats = dset.variables['XLAT_M'][0]
            lons = dset.variables['XLONG_M'][0]
            hgt = dset.variables['HGT_M'][0]
        xhgt, yhgt = m(lons, lats)
        hgt = gaussian_filter(hgt, 1)
        return xhgt, yhgt, hgt
    xhgt, yhgt, hgt = get_terrain()
    CS = plt.contour(xhgt, yhgt, hgt, np.arange(500,5000,500), linewidths=1, colors='k', alpha = 0.5)

    plt.colorbar(cs)
    ax.set_title(var)

    plt.savefig('./surf_figs/' +var)
    plt.close()
