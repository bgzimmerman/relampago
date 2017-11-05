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

os.chdir('/home/disk/meso-home/bzim/relampago')

cormenlat = (-31.42,-32.89)
cormenlon = (-64.18,-68.84)
exportdir = '/home/disk/user_www/bzim/relampago'
variables_gfs = {'Precipitable_water_entire_atmosphere_single_layer' : ('Precipitable Water','mm','pwat',(0,35),'ncar_moist'),
             'Precipitation_rate_surface_Mixed_intervals_Average' : ('Total Precipitation Rate','mm/hr','precip',(0,10), 'ncar_precip'),
'Apparent_temperature_height_above_ground' : ('2m Temperature', 'degK', '2mtemp', (260, 312), 'ncar_temp'),
'Relative_humidity_height_above_ground' : ('2m Relative Humidity', '%', '2mRH', (0, 100), 'MPL_RdBu'),
'Convective_available_potential_energy_surface': ('Surface CAPE', 'J/Kg', 'cape', (0, 4000), 'ncar_ua')
             }

# Spread variables need to have a range,
# Mean variables need to have a tuple
variables_gefs = {'Convective_available_potential_energy_pressure_difference_layer_unweightedMean' :
                  ('CAPE','J/kg','cape',(0,4000,100),'ncar_ua'),
             'Convective_available_potential_energy_pressure_difference_layer_stdDev' :
                  ('CAPE','J/kg','cape',np.arange(0,600), 'ncar_precip'),
             'Total_precipitation_surface_6_Hour_Accumulation_unweightedMean' :
                  ('6-hr Precipitation','mm','precip06',(0,30,1), 'ncar_precip'),
             'Total_precipitation_surface_6_Hour_Accumulation_stdDev' :
                  ('6-hr Precipitation','mm','precip06',np.arange(15), 'ncar_precip'),
             'Pressure_reduced_to_MSL_msl_unweightedMean' :
                  ('Sea-Level Pressure','hPa','slp',(940,1060, 4), 'ncar_precip'),
             'Pressure_reduced_to_MSL_msl_stdDev' :
                  ('Sea-Level Pressure Sprd',
                   'hPa','slp',np.arange(20),'WhBlGrYeRe'),
             'Geopotential_height_isobaric_unweightedMean' :
                  ('500 hPa Height','m','hgt500',(5400,6060, 60), 'ncar_precip'),
             'Geopotential_height_isobaric_stdDev' :
                  ('500 hPa Height Spread','m','hgt500',np.arange(30), 'ncar_precip'),
             'Relative_humidity_height_above_ground_unweightedMean' :
                    ('2m Relative Humidity', '%', '2mRH', (0, 100, 5), 'ncar_wv'),
            'Relative_humidity_height_above_ground_stdDev' :
                    ('2m Relative Humidity Spread', '%', '2mRH', np.arange(20), 'ncar_wv'),
             'Temperature_height_above_ground_unweightedMean' :
                    ('2m Temperature', 'degK', '2mtemp', (-15, 40, 1), 'ncar_temp'),
            'Temperature_height_above_ground_stdDev' :
                    ('2m Temperature Spread', 'degK', '2mtemp', np.arange(5), 'ncar_temp'),


                 }





def get_gfs(type='hires', bounds=(-72.0, -54.5, -37.0, -26.25)):
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
    if type=='hires':
        query.variables(*((variables_gfs.keys())))
    elif type=='ensemble':
        query.variables(*((variables_gefs.keys())))
    data = ncss.get_data(query)
    return data





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


# In[146]:

def make_plot(pltdat, export=False):
    gs = gridspec.GridSpec(2,2,height_ratios=(25,1))
    plt.figure(figsize=(10,8))
    plt.subplot(gs[0,:])

    if pltdat['pltcode'] == 'precip':
        pltdat['field'] = pltdat['field'] * 3600.0
    if pltdat['pltcode'] in ['slp','hgt500', 'hgt250', 'hgt925']:
        print(pltdat['field'].min(), pltdat['field'].max())
        print(pltdat['contourfield'].min(), pltdat['contourfield'].max())
        lw=2
    else:
        lw=2
    if pltdat['pltcode'] == '2mtemp':
        key = 'contourfield' if pltdat['modelname'] == 'GEFS' else 'field'
        pltdat[key] = pltdat[key]- 273.15

    if 'contourfield' in pltdat.keys():
        spreadfield = gaussian_filter(pltdat['contourfield'],1)
        lines = plt.contour(x,y,spreadfield, np.arange(pltdat['contourlevs'][0], pltdat['contourlevs'][1], pltdat['contourlevs'][2]), linewidths=lw, colors = 'k')
        clabels = plt.clabel(lines, fmt='%2.0f')

    if pltdat['modelname']=='GEFS':
        theplot = plt.pcolormesh(x,y,pltdat['field'], vmin=pltdat['range'][0], vmax=pltdat['range'][-1], cmap=color_map('WhBlGrYeRe'))
    elif pltdat['modelname']=='GFS':
        theplot = plt.pcolormesh(x,y,pltdat['field'], vmin=pltdat['range'][0], vmax=pltdat['range'][-1], cmap=pltdat['cmap'])
    # OLD theplot = plt.pcolormesh(x,y,pltdat['field'], vmin=pltdat['range'][0], vmax=pltdat['range'][1], cmap=pltdat['cmap'])
    # Highlight key sites
    plt.scatter(xpt, ypt, s=100, marker='s', facecolor='Firebrick', edgecolor='white')

    # Contour terrain
    if pltdat['plot terrain']:
        plt.contour(xhgt, yhgt, hgt, np.arange(500,4000,500), linewidths=1, colors='0.4')
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
            os.system('mkdir {:s}'.format(outpath))

        plt.savefig(filename, bbox_inches='tight')
        os.system('mv {:s} {:s}'.format(filename, outpath))

        if not os.path.exists(outpath2):
            os.system('mkdir {:s}'.format(outpath2))
        print(filename2)
        plt.savefig(filename2, bbox_inches='tight')
        os.system('mv {:s} {:s}'.format(filename2, outpath2))
        os.system('convert {:s}/{:s} {:s}/{:s}'.format(outpath2, filename2, outpath2, filename2[:-3]+'gif'))
        os.system('rm {:s}/{:s}'.format(outpath2, filename2))
        plt.close()
    else:
        plt.show()

if __name__ == '__main__':

    # now for the ensemble
    # Also get lat/lon bounds
    ensbounds = (-110.0, -46.0, -50.0, -11.0)
    m = Basemap(projection='merc', llcrnrlon=ensbounds[0],
                urcrnrlon=ensbounds[1], llcrnrlat=ensbounds[2],
                urcrnrlat=ensbounds[3], resolution='l', area_thresh=5000)

    #xhgt, yhgt, hgt = get_terrain()

    data = get_gfs(type='ensemble', bounds=ensbounds)
    glon, glat = np.meshgrid(data.variables['lon'][:]-360, data.variables['lat'][:])
    x,y = m(glon, glat)
    xpt, ypt = m(cormenlon, cormenlat)

    if 'time1' in data.variables.keys():
        timevar = 'time1'
    elif 'time2' in data.variables.keys():
        timevar = 'time2'
    elif 'time8' in data.variables.keys():
        timevar = 'time8'
    else:
        timevar = 'time'
    times = num2date(data.variables[timevar][:], units = data.variables[timevar].units)
    for n,time in enumerate(times):
        #if n !=0:
        #    continue
        print(time)
        flead = data.variables[timevar][n]
        init = time - timedelta(hours=flead)
        print(init)
        #print(init)
        for fullvar, varcodes in variables_gefs.items():
            if fullvar.endswith('stdDev'):
                continue
            varpts = fullvar.split('_')
            varpts[-1] = 'stdDev'
            stdevvar = '_'.join(varpts)

            print('   ', varcodes[2])

            if varcodes[2] in ['slp']:
                # Contour the mean
                meanlevs = varcodes[3]
                varcodes = variables_gefs[stdevvar]
                plotdata = {
                'field' : np.squeeze(data.variables[stdevvar][n]/100.0),
                'contourfield' : np.squeeze(data.variables[fullvar][n]/100.0),
                'contourlevs' : meanlevs,
                'varname' : varcodes[0],
                'varunit' : varcodes[1],
                'modelname' : 'GEFS',
                'cmap' : color_map(varcodes[4]),
                'valid' : time,
                'init' : init,
                'flead' : int(flead),
                'pltcode' : varcodes[2],
                'range' : varcodes[3],
                'plot terrain' : False,
                }
                make_plot(plotdata, export=True)

            elif varcodes[2].startswith('hgt'):
                #print(data.variables[fullvar])
                #import pdb; pdb.set_trace()
                for level in ['500']:
                    try:
                        isolevs = data.variables['isobaric1'][:]
                    except:
                        isolevs = data.variables['isobaric2'][:]
                    thislev = int(level)*100
                    hdex = np.argmax(isolevs >= thislev)
                    #print(isolevs[hdex])
                    #exit()
                    # Contour the mean
                    meanlevs = varcodes[3]
                    varcodes = variables_gefs[stdevvar]
                    plotdata = {
                    'field' : np.squeeze(data.variables[stdevvar][n, hdex]),
                    'contourfield' : np.squeeze(data.variables[fullvar][n, hdex]),
                    'contourlevs' : meanlevs,
                    'varname' : varcodes[0],
                    'varunit' : varcodes[1],
                    'modelname' : 'GEFS',
                    'cmap' : color_map(varcodes[4]),
                    'valid' : time,
                    'init' : init,
                    'flead' : int(flead),
                    'pltcode' : varcodes[2],
                    'range' : varcodes[3],
                    'plot terrain' : False,
                    }
                    """plotdata = {
                    'field' : np.squeeze(data.variables[stdevvar][n,hdex]),
                    'contourfield' : np.squeeze(data.variables[fullvar][n,hdex]),
                    'contourlevs' : meanlevs,
                    'varname' : varcodes[0],
                    'varunit' : varcodes[1],
                    'modelname' : 'GEFS',
                    'cmap' : color_map(varcodes[4]),
                    'valid' : time,
                    'init' : init,
                    'flead' : int(flead),
                    'pltcode' : 'hgt' + level,
                    'range' : datarange,
                    'plot terrain' : False,
                    }"""
                    make_plot(plotdata, export=True)

            else:
                # Contour the mean
                meanlevs = varcodes[3]
                varcodes = variables_gefs[stdevvar]
                plotdata = {
                'field' : np.squeeze(data.variables[stdevvar][n]),
                'contourfield' : np.squeeze(data.variables[fullvar][n]),
                'contourlevs' : meanlevs,
                'varname' : varcodes[0],
                'varunit' : varcodes[1],
                'modelname' : 'GEFS',
                'cmap' : color_map(varcodes[4]),
                'valid' : time,
                'init' : init,
                'flead' : int(flead),
                'pltcode' : varcodes[2],
                'range' : varcodes[3],
                'plot terrain' : False,
                }
                """
                # Contour the stdev
                plotdata = {
                'field' : np.squeeze(data.variables[fullvar][n]),
                'contourfield' : np.squeeze(data.variables[stdevvar][n]),
                'contourlevs' : variables_gefs[stdevvar][3],
                'varname' : varcodes[0],
                'varunit' : varcodes[1],
                'modelname' : 'GEFS',
                'cmap' : color_map(varcodes[4]),
                'valid' : time,
                'init' : init,
                'flead' : int(flead),
                'pltcode' : varcodes[2],
                'range' : varcodes[3],
                'plot terrain' : False,
                }
                """
                make_plot(plotdata, export=True)


    # First, hires
    bounds=(-72.0, -54.5, -37.0, -26.25)
    m = Basemap(projection='merc', llcrnrlon=bounds[0], urcrnrlon=bounds[1], llcrnrlat=bounds[2], urcrnrlat=bounds[3], resolution='l', area_thresh=5000)

    xhgt, yhgt, hgt = get_terrain()
    data = get_gfs(type='hires', bounds=bounds)
    glon, glat = np.meshgrid(data.variables['lon'][:]-360, data.variables['lat'][:])
    x,y = m(glon, glat)
    xpt, ypt = m(cormenlon, cormenlat)

    if 'time1' in data.variables.keys():
        timevar = 'time1'
    elif 'time2' in data.variables.keys():
        timevar = 'time2'
    elif 'time8' in data.variables.keys():
        timevar = 'time8'
    else:
        timevar = 'time'

    times = num2date(data.variables[timevar][:], units = data.variables[timevar].units)
    for n,time in enumerate(times):
        #if n !=0:
        #    continue
        print(time)
        flead = data.variables[timevar][n]
        init = time - timedelta(hours=flead)
        print(init)
        #print(init)
        for fullvar, varcodes in variables_gfs.items():
            print('   ', varcodes[2])
            plotdata = {
            'field' : data.variables[fullvar][n].squeeze(),
            'varname' : varcodes[0],
            'varunit' : varcodes[1],
            'modelname' : 'GFS',
            'cmap' : color_map(varcodes[4]),
            'valid' : time,
            'init' : init,
            'flead' : int(flead),
            'pltcode' : varcodes[2],
            'range' : varcodes[3],
            'plot terrain' : True,
            }
            make_plot(plotdata, export=True)
