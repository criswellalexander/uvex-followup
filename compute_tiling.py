#!/usr/bin/env python
# coding: utf-8

from astropy import units as u
from ligo.skymap.tool import ArgumentParser, FileType
from dorado.scheduling import mission as _mission
from dorado.scheduling.units import equivalencies

from astropy.coordinates import ICRS
from astropy_healpix import HEALPix
from astropy.io import fits
from astropy.time import Time
from astropy.table import QTable
from ligo.skymap.io import read_sky_map
from ligo.skymap.bayestar import rasterize
from ligo.skymap import plot
from ligo.skymap.postprocess import find_greedy_credible_levels
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from tqdm import tqdm
from matplotlib.ticker import MaxNLocator
import sys, os, configparser
import argparse
import pandas as pd

def compute_tiling(allsky_sched,fitsloc,schedloc,outdir,duration):
    '''
    Function to compute the number of UVEX pointings needed to tile a GW localization region.
    
    Arguments
    -----------------------------
    allsky_sched (str)    : /path/to/allsky_scheduled.txt
    fitsloc (str)         : /path/to/directory/with/skymaps/
    schedloc (str)        : /path/to/directory/with/schedules/ (as computed by the uvex scheduler)
    outdir (str)          : /path/to/save/directory/
    duration (str)        : observing duration in hours, written as (e.g.) '3 hr'
    '''
    ## manually define some arguments
    duration = u.Quantity(duration)
    time_step = u.Quantity('1 min')
#     mission = getattr(_mission, 'uvex')
    mission = _mission.uvex()
    nside = 64
    delay = 0

    healpix = HEALPix(nside, order='nested', frame=ICRS())

    ## get event numbers
    event_df = pd.read_csv(allsky_sched,delimiter=' ')#,skiprows=1)

    event_list = event_df['event_id'].tolist()

    texp_list = event_df['t_exp (ks)'].tolist()

    ## get paths to data file directories
    fitslist = list(map(lambda ev_num: str(ev_num)+'.fits',event_list))
    schedlist = list(map(lambda ev_num: str(ev_num)+'.ecsv',event_list))


    rows = []

    ## loop over .fits and .ecsv files
    for event, fitsfile, schedule, texp in zip(event_list,fitslist,schedlist,texp_list):

        # Read multi-order sky map and rasterize to working resolution
        skymap = read_sky_map(fitsloc+fitsfile, moc=True)['UNIQ', 'PROBDENSITY']
        skymap = rasterize(skymap, healpix.level)['PROB']
        # check to see if file is empty before loading because there are empty files for some reason
        if os.stat(schedloc+schedule).st_size == 0:
            print("Error: empty schedule for event",event)
            continue    
        schedule = QTable.read(schedloc+schedule, format='ascii.ecsv')

        indices = np.asarray([], dtype=np.intp)
        tiles_to_99pct = None
        row_count = 0
        reached_99 = False

        for row in schedule:
            row_count += 1
            new_indices = mission.fov.footprint_healpix(
                healpix, row['center'], row['roll'])
            indices = np.unique(np.concatenate((indices, new_indices)))
            if (100*skymap[indices].sum() > 99) and (reached_99==False):
                tiles_to_99pct = row_count
                reached_99 = True
        tiles_total = row_count
        total_prob = 100*skymap[indices].sum()
        rows.append([event,total_prob,texp,tiles_to_99pct,tiles_total])

        print("Percent coverage is",total_prob,"% for event",event)

    ## Format and save
    events_cov = pd.DataFrame(rows,columns=['event_id','percent_coverage','texp_sched (ks)','tiles_to_99pct','tiles_total'])
    # events_cov['texp_sched (s)'] = texp_list
#     outname = allsky_sched.split('_')[-1].split('.')[0]
    events_cov.to_csv(outdir+'/allsky_coverage.txt',index=False,sep=' ')
    
    return


if __name__ == '__main__':
    
    ## set up argparser
    parser = argparse.ArgumentParser(description="Given observing schedules, compute tiling and statistics.")
    parser.add_argument('params', type=str, help='/path/to/params_file.ini')
    
    args = parser.parse_args()
    
    ## set up configparser
    config = configparser.ConfigParser()
    config.read(args.params)
    
    ## get info from params file
    obs_scenario_dir = config.get("params","obs_scenario")
    out_dir          = config.get("params","save_directory")
    duration         = config.get("params","tiling_time")
    
    allsky_sched = out_dir+'/allsky_sched_full.txt'
    fitsloc = obs_scenario_dir+'/allsky/'
    schedloc = out_dir+'/schedules/'
    
    ## run the script
    compute_tiling(allsky_sched,fitsloc,schedloc,out_dir,duration)