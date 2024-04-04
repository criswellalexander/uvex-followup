#!/usr/bin/env python
# coding: utf-8

# Not sure how this copyright works...
# Adapted by Alexander Criswell from animate.py in the NASA Dorado code for the purpose of calculating
# tiling considerations for the UVEX mission proposal
#
# Copyright © 2021 United States Government as represented by the Administrator
# of the National Aeronautics and Space Administration. No copyright is claimed
# in the United States under Title 17, U.S. Code. All Other Rights Reserved.
#
# SPDX-License-Identifier: NASA-1.3
#
"""Compute the number of UVEX pointings to cover a GW localization region of <100 sq deg"""

'''
Usage: python UVEX_Tiling_Calculator.py [/path/to/allsky_texp.txt] [/path/to/fits/directory/] [/path/to/schedules/directory/] [/path/to/output/directory/]
'''

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
import sys
import os
import pandas as pd

allsky_file = str(sys.argv[1])
fitsloc = str(sys.argv[2])
schedloc = str(sys.argv[3])
outdir = str(sys.argv[4])

## manually define some arguments
duration = u.Quantity('3 hr')
time_step = u.Quantity('1 min')
mission = getattr(_mission, 'uvex')
nside = 64
delay = 0

healpix = HEALPix(nside, order='nested', frame=ICRS())

## get event numbers
event_df = pd.read_csv(allsky_file,delimiter=' ')#,skiprows=1)

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
outname = allsky_file.split('_')[-1].split('.')[0]
events_cov.to_csv(outdir+'/allsky_coverage_'+outname+'.txt',index=False,sep=' ')
