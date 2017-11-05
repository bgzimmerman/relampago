from datetime import datetime, timedelta
import requests
import os
EV = os.environ

# This script will download text files from Paola's website

current_time = datetime.now()

def download_ARsurf(sample_time):

    base_url = 'http://cimaps.cima.fcen.uba.ar/relampago/estaciones'
    full_url = '{:s}/{:%Y/%m/%d/%Y%m%d%H}.csv'.format(base_url, sample_time)
    outpath = '{:s}/relampago/AR_surf_data/{:%Y%m%d%H}.csv'.format(
                EV['HOME'], sample_time)

    response = requests.get(full_url)
    with open(outpath, 'w') as f:
        f.write(response.content)
    print '{:%Y%m%d%H} downloaded'.format(sample_time)
    return

def get_all_files():
    # This function downloads all the files since the start of
    # data collection

    # Set up the starting date and the time increment
    startdate = datetime(2017, 10, 01, 0)
    hour = timedelta(1./24)

    # Increment startdate by one hour at a time and get the data
    date = startdate
    while date <= datetime.now():
        download_ARsurf(date)
        date += hour

if __name__ == '__main__':
    download_ARsurf(datetime.now())
    #get_all_files()
