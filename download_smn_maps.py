#!/home/disk/aleutians/nobackup/lmadaus/anaconda/bin/python
from __future__ import print_function, division
from datetime import datetime, timedelta
import urllib2
import urllib
#import urllib.request
import os

os.chdir('/home/disk/aleutians/lmadaus/relampago')

def download_images_smn():
    imgdir = '/home/disk/user_www/lmadaus/relampago/smneta'
    webpath = 'http://www.smn.gov.ar/dpd/cartas'

    filetemp = {'slp' : ('slp_{:d}.gif', [1,9,17,25,33]),
                'precip24' : ('ppp24dia00_{:d}.gif', [24,48,72,96]),
                't2rh' : ('th_{:d}.gif', [5,7,13,15,21,23]),
                'u10rdp' : ('rdp_color_{:d}.gif', [1,3,5,6,9,11,13,15,17]),
                'precip03' : ('pppv_{:d}.gif', [2,3,4,5,6,7,8,9])}
    for imgname, imginfo in filetemp.items():
        print(imgname)
        imgtemp = imginfo[0]

        for imgcount, imgnum in enumerate(imginfo[1]):
            fullpath = '/'.join((webpath, imgtemp.format(imgnum)))
            # Check the past updated time if this is the first image
            if imgcount == 0:
                # Python 2 interface
                conn = urllib2.urlopen(fullpath, timeout=30)
                last_modified = conn.info().getdate('last-modified')
                imgtime = datetime(*last_modified[:-2])
                # Back out the model init time
                while imgtime.hour not in [0,12]:
                    imgtime -= timedelta(hours=1)
                modelinit = datetime(imgtime.year, imgtime.month, imgtime.day, imgtime.hour)
                # Python 3 interface
                #last_modified = conn.header['last-modified']
                print('    ', modelinit)

                # check to see if this init time exists in web dir
                outpath ='/'.join((imgdir, modelinit.strftime('%Y%m%d%H')))
                if not os.path.exists(outpath):
                    os.system('mkdir {:s}'.format(outpath))
            
            # figure out the valid date
            valid = modelinit + timedelta(hours=imgnum*3)
            # Now actually download the file
            fulloutpath = '/'.join((outpath, '{:s}_{:%Y%m%d%H}_f{:03d}.gif'.format(imgname, valid, imgnum*3)))
            urllib.urlretrieve(fullpath, fulloutpath)

def download_images_dmc():
    imgdir = '/home/disk/user_www/lmadaus/relampago/dmcwrf'
    webpath = 'http://archivos.meteochile.gob.cl/portaldmc/meteochile/wrf/chile'

    filetemp = {'slp' : ('slp_espesor1000-500_D1_{:03d}.png', range(1,42)),
                'precip03' : ('slp_espesor1000-500_prec_D1_{:03d}.png', range(1,41)),
                'u10rainc' : ('Viento10m_PCONV_wrf_D1_{:03d}.png', range(1,41)),
                'capepwat' : ('cape_aguaprec_D1_{:03d}.png', range(1,42)),
                'vort500' : ('vort_geo_tmp_500_D1_{:03d}.png', range(1,42)),
               }
    for imgname, imginfo in filetemp.items():
        print(imgname)
        imgtemp = imginfo[0]

        for imgcount, imgnum in enumerate(imginfo[1]):
            fullpath = '/'.join((webpath, imgtemp.format(imgnum)))
            #print(fullpath)
            # Check the past updated time if this is the first image
            if imgcount == 0:
                # Python 2 interface
                conn = urllib2.urlopen(fullpath, timeout=30)
                last_modified = conn.info().getdate('last-modified')
                imgtime = datetime(*last_modified[:-2])
                # Back out the model init time
                while imgtime.hour not in [0,12]:
                    imgtime -= timedelta(hours=1)
                modelinit = datetime(imgtime.year, imgtime.month, imgtime.day, imgtime.hour)
                # Python 3 interface
                #last_modified = conn.header['last-modified']
                print('    ', modelinit)

                # check to see if this init time exists in web dir
                outpath ='/'.join((imgdir, modelinit.strftime('%Y%m%d%H')))
                if not os.path.exists(outpath):
                    os.system('mkdir {:s}'.format(outpath))
            
            # figure out the valid date
            if imgname in ['precip03','u10rainc']:
                flead = (imgnum+1) * 3
            else:
                flead = imgnum * 3
            valid = modelinit + timedelta(hours=flead)
            # Now actually download the file
            fulloutpath = '/'.join((outpath, '{:s}_{:%Y%m%d%H}_f{:03d}.png'.format(imgname, valid, flead)))
            urllib.urlretrieve(fullpath, fulloutpath)

def download_images_inmet():
    imgdir = '/home/disk/user_www/lmadaus/relampago/inmetcosmo'
    webpath = 'http://www.inmet.gov.br/projetos/cga/capre/cosmo7/AS'

    filetemp = {
        'precip24' : ('prec24h/web_AS_precip24h_{:%Y%m%d%H}00_+{:d}.png', range(24,175,24)),
        'precip03' : ('prec3hPressV10m/web_AS_Prec3hPressV10m_{:%Y%m%d%H}00_+{:d}.png', range(3,175,3)),
        'capek' : ('Kcape/web_AS_Kcape_{:%Y%m%d%H}00_+{:d}.png', range(0,175,3)),
        'h500' : ('NuvemMediaGeoTemp500/web_AS_NuvemMediaGeoTemp500_{:%Y%m%d%H}00_+{:d}.png', range(0,175,3)),
        'pwat850' : ('aguaprecT850V850/web_AS_aguaprecT850V850_{:%Y%m%d%H}00_+{:d}.png', range(0,175,3)),
               }

    # Here we have to guess on modelinit
    now = datetime.utcnow()
    while now.hour not in [0,12]:
        now -= timedelta(hours=1)
    modelinit = datetime(now.year, now.month, now.day, now.hour)
    modelinit -= timedelta(hours=6)
    for imgname, imginfo in filetemp.items():
        print(imgname)
        imgtemp = imginfo[0]

        for imgcount, imgnum in enumerate(imginfo[1]):
            fullpath = '/'.join((webpath, imgtemp.format(modelinit, imgnum)))
            #print(fullpath)
            # Check the past updated time if this is the first image
            if imgcount == 0:
                # check to see if this init time exists in web dir
                outpath ='/'.join((imgdir, modelinit.strftime('%Y%m%d%H')))
                if not os.path.exists(outpath):
                    os.system('mkdir {:s}'.format(outpath))
            
            flead = imgnum 
            valid = modelinit + timedelta(hours=flead)
            # Now actually download the file
            fulloutpath = '/'.join((outpath, '{:s}_{:%Y%m%d%H}_f{:03d}.png'.format(imgname, valid, flead)))
            #print(fullpath)
            urllib.urlretrieve(fullpath, fulloutpath)

if __name__ == '__main__':
   #download_images_smn()
   #download_images_dmc()
   download_images_inmet()



