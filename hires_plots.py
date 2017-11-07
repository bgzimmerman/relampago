from utils import *

os.chdir(EV['HOME'] + '/relampago')
exportdir = EV['HOME'] + '/relampagofigs'
# Build the set of variables to be explored by adding dictionary entries
exportdir = '/home/disk/user_www/bzim/relampago'
# Unweighted Means have contour_intervals, while stdDev has spread_range
# Potentially put these in later?
variables_gefs = {}
variables_gefs['Geopotential_height_isobaric_unweightedMean'] = {
                'short_name' : '{:.0f} hPa Height',
                'units' : 'm',
                'pltcode' : 'hgt{:.0f}',
                'contour_info' : {
                    '500' : ((4000, 7000, 60), np.arange(100)),
                    '925' : ((0, 4000, 30), np.arange(100)),
                    '250' : ((8000, 13000, 120), np.arange(120))
                    }
                'cmap' : 'ncar_precip'}

variables_gefs['Geopotential_height_isobaric_stdDev'] = {
                'short_name' : '{:.0f} hPa Height Spread',
                'units' : 'm',
                'pltcode' : 'hgt{:.0f}',
                'cmap' : 'ncar_precip'}

variables_gefs['Temperature_isobaric_unweightedMean'] = {
                'short_name' : '{:.0f} hPa Temperature',
                'units' : 'degC',
                'pltcode' : 'temp{:.0f}',
                'contour_info' : {
                    '500' : ((0, 330, 2), np.arange(10)),
                    '925' : ((0, 330, 2), np.arange(10)),
                    '250' : ((0, 330, 2), np.arange(10))
                    }
                'cmap' : 'ncar_temp'}

variables_gefs['Temperature_isobaric_stdDev'] = {
                'short_name' : '{:.0f} hPa Temperature Spread',
                'units' : 'degC',
                'pltcode' : 'temp{:.0f}',
                'cmap' : 'ncar_temp'}

ensbounds = (-110.0, -46.0, -50.0, -11.0)
m = Basemap(projection='merc', llcrnrlon=ensbounds[0],
            urcrnrlon=ensbounds[1], llcrnrlat=ensbounds[2],
            urcrnrlat=ensbounds[3], resolution='l', area_thresh=5000)

#xhgt, yhgt, hgt = get_terrain()

data = get_gfs(variables_gefs, type='ensemble', bounds=ensbounds)
glon, glat = np.meshgrid(data.variables['lon'][:]-360, data.variables['lat'][:])

for fullvar, varcodes in variables_gefs.items():
    if fullvar.endswith('stdDev'):
        # This is because  the stdDev and unweightedMean need to be
        # put on the same plot, and we're trying to iterate through the items.
        # all that's really different though is the short name and the contour_intervals
        # intervals and spread_range... May be a better way?
        continue
    # Make stdevvar from fullvar
    varpts = fullvar.split('_')
    varpts[-1] = 'stdDev'
    stdevvar = '_'.join(varpts)

    print(varcodes['pltcode'])

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
        print(time)
        flead = data.variables[timevar][n]
        init = time - timedelta(hours=flead)
        print(init)
        isolevs = data.variables['isobaric2'][:]
        for level in [92500., 50000., 25000.]:
            hdex = isolevs == level
            plotdata = {
            'field' : np.squeeze(data.variables[stdevvar][n, hdex]),
            'contourfield' : np.squeeze(data.variables[fullvar][n, hdex]),
            'contourlevs' : varcodes['contour_info']['{:.0f}'.format(level/100)][0],
            'glon' : glon,
            'glat' : glat,
            'varname' : varcodes['short_name'].format(level/100),
            'varunit' : varcodes['units'],
            'modelname' : 'GEFS',
            'cmap' : color_map(varcodes['cmap']),
            'valid' : time,
            'init' : init,
            'flead' : int(flead),
            'pltcode' : varcodes['pltcode'].format(level/100),
            'range' : varcodes['contour_info']['{:.0f}'.format(level/100)][1],
            'plot terrain' : False,
            }
            make_plot(plotdata, m, export=True)
