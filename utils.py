#!usr/bin/python

from __future__ import print_function, division
from datetime import datetime
import numpy as np
from datetime import datetime, timedelta
from netCDF4 import Dataset, num2date
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
#plt.style.use('ncar_wrf_plots')
plt.style.use('ggplot')
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from color_maker.color_maker import color_map
from matplotlib import cm
import os
from scipy.ndimage.filters import gaussian_filter
from os import environ as EV
exportdir = EV['HOME'] + '/relampagofigs'


def get_gfs(variables, type='hires', bounds=(-72.0, -54.5, -37.0, -26.25)):
    if type == 'hires':
        best_gfs = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/latest.xml')
    elif type=='ensemble':
        best_gfs = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/NCEP/GEFS/Global_1p0deg_Ensemble/derived/latest.xml')

    best_ds = list(best_gfs.datasets.values())[0]
    ncss = NCSS(best_ds.access_urls['NetcdfSubset'])
    query = ncss.query()
    # West east south north
    query.lonlat_box(*bounds)
    query.accept('netcdf4')
    query.time_range(datetime.utcnow(), datetime.utcnow()+timedelta(days=5))
    query.variables(*((variables.keys())))
    data = ncss.get_data(query)
    return data

def get_terrain(terrfile=EV['HOME'] + 'relampago/geo_em.d01.nc'):
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

def make_plot(pltdat, m, export=False):
    lw = 2
    gs = gridspec.GridSpec(2,2,height_ratios=(25,1))
    plt.figure(figsize=(10,8))
    plt.subplot(gs[0,:])

    if pltdat['pltcode'] == 'precip':
        pltdat['field'] = pltdat['field'] * 3600.0

    """
    if pltdat['pltcode'] in ['slp','hgt500', 'hgt250', 'hgt925']:
        print(pltdat['field'].min(), pltdat['field'].max())
        print(pltdat['contourfield'].min(), pltdat['contourfield'].max())
        lw=2
    else:
        lw=2
    """

    if pltdat['pltcode'] == '2mtemp':
        key = 'contourfield' if pltdat['modelname'] == 'GEFS' else 'field'
        pltdat[key] = pltdat[key]- 273.15

    if 'contourfield' in pltdat.keys():
        spreadfield = gaussian_filter(pltdat['contourfield'],1)
        lines = m.contour(pltdat['glon'],pltdat['glat'],spreadfield, np.arange(pltdat['contourlevs'][0], pltdat['contourlevs'][1], pltdat['contourlevs'][2]), linewidths=lw, colors = 'k', latlon=True)
        clabels = plt.clabel(lines, fmt='%2.0f')

    if pltdat['modelname']=='GEFS':
        theplot = m.pcolormesh(pltdat['glon'],pltdat['glat'],pltdat['field'], vmin=pltdat['range'][0], vmax=pltdat['range'][-1], cmap=color_map('WhBlGrYeRe'), latlon=True)
    elif pltdat['modelname']=='GFS':
        theplot = m.pcolormesh(pltdat['glon'],pltdat['glat'],pltdat['field'], vmin=pltdat['range'][0], vmax=pltdat['range'][-1], cmap=pltdat['cmap'], latlon=True)
    # OLD theplot = plt.pcolormesh(x,y,pltdat['field'], vmin=pltdat['range'][0], vmax=pltdat['range'][1], cmap=pltdat['cmap'])
    # Highlight key sites
    #plt.scatter(xpt, ypt, s=100, marker='s', facecolor='Firebrick', edgecolor='white')

    # Contour terrain
    if pltdat['plot terrain']:
        m.contour(xhgt, yhgt, hgt, np.arange(500,4000,500), linewidths=1, colors='0.4')
    m.drawcoastlines(linewidth=2)
    m.drawcountries(linewidth=2)
    m.drawstates(linewidth=1, color='0.8')
    parallels = np.arange(-60,-5, 5)
    meridians = np.arange(-120,-30, 5)
    m.drawparallels(parallels,labels=[False,True,True,False])
    m.drawmeridians(meridians,labels=[True,False,True,False])

    tax = plt.subplot(gs[1,0])
    plt.text(0,0.75, '{modelname:s} {varname:s}'.format(**pltdat), fontsize=16, ha='left', va='center', transform=tax.transAxes)
    plt.text(0,0.00, 'Valid: {valid:%d %b %H00 UTC} | Init: {init:%d %b %H00 UTC}'.format(**pltdat), fontsize=12, ha='left', va='center', transform=tax.transAxes)
    plt.text(1.0, 0.75, 'F{flead:03d}'.format(**pltdat), fontsize=16, ha='right', va='center', transform=tax.transAxes)
    #plt.text(0,-0.75, 'Init: {init: %d %b %H00 UTC}'.format(**pltdat), fontsize=18, ha='left', va='center', transform=tax.transAxes)
    tax.axis('off')

    cax = plt.subplot(gs[1,1])
    cbar = plt.colorbar(theplot, cax=cax, orientation='horizontal')
    cbar.set_label(label='{varname:s} [{varunit:s}]'.format(**pltdat), fontsize=12)
    cax.tick_params(labelsize=18)
    plt.tight_layout()

    if export:
        outpath = '{:s}/{:s}/{:%Y%m%d%H}'.format(exportdir,pltdat['modelname'].lower(), pltdat['init'])
        filename = '{pltcode:s}_{valid:%Y%m%d%H}_f{flead:03d}.png'.format(**pltdat)

        outpath2= EV['HOME']+'/relampago/ftp_gifs'
        filename2= 'model.{modelname:s}.{init:%Y%m%d%H}.{flead:03d}_{pltcode:s}.png'.format(**pltdat)

        if not os.path.exists(outpath):
            os.system('mkdir -p {:s}'.format(outpath))

        plt.savefig(filename, bbox_inches='tight')
        os.system('mv {:s} {:s}'.format(filename, outpath))

        if not os.path.exists(outpath2):
            os.system('mkdir -p {:s}'.format(outpath2))
        print(filename2)
        plt.savefig(filename2, bbox_inches='tight')
        os.system('mv {:s} {:s}'.format(filename2, outpath2))
        os.system('convert {:s}/{:s} {:s}/{:s}'.format(outpath2, filename2, outpath2, filename2[:-3]+'gif'))
        os.system('rm {:s}/{:s}'.format(outpath2, filename2))
        plt.close()
    else:
        plt.show()
