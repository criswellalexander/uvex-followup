#!/usr/bin/env python
# coding: utf-8

# This script contains some code to:
# 1. Determine the maximum exposure time across a GW sky localization map for all events with a <100 square degree localization.
# 2. Discard all events with t_exp_max > 10 ks
# 3. Save a .csv with the event number and t_exp_mean
# It assumes a fiducial absolute KNe AB magnitude of -12.1 and uses the mean GW distance estimate (DISTMEAN) to estimate the resulting apparent magnitude.

'''
Usage: texp-by-sky-loc.py [/path/to/directory/with/skymaps/] [/path/to/output/directory/] [/path/to/allsky.txt]

Arguments:
    [1] path to location of LIGO GW localization skymap files.
    [2] output directory. Needs to already exist.
    [3] CSV describing GW simulations. At minimum needs to have 'simulation_id', 'area_(90)', and 'distmean' columns
    [4] band ('nuv' or 'fuv')
    [5] M_AB in that band
'''

from functools import partial

from astropy_healpix import HEALPix
from astropy.coordinates import ICRS
from astropy import units as u
from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
import healpy as hp
import pandas as pd

from synphot import SourceSpectrum, Observation, SpectralElement
from synphot.models import ConstFlux1D, BlackBodyNorm1D

import astropy.units as u
import numpy as np

import uvex.sensitivity
from uvex.sensitivity import zodi, galactic

from importlib import reload
reload(uvex.sensitivity)
reload(zodi)

## vim is the absolute worst and I can't copy/paste this lower without making a mess
## so enjoy this random nonsensical out-of-order import
import sys
band = sys.argv[4]

nuv_band = uvex.sensitivity.filters.nuv_bandpass()
fuv_band = uvex.sensitivity.filters.fuv_bandpass()
if band=='nuv':
    wave = np.arange(1300, 10000)*u.AA
elif band=='fuv':
    wave = np.arange(1000, 10000)*u.AA
area = uvex.sensitivity.config.AREA

from astropy.coordinates import SkyCoord
from astropy.time import Time# Pick your favorite date
obstime = Time('2021-02-18 09:00:00')

from astropy import units as u
from ligo.skymap.tool import ArgumentParser, FileType
from dorado.scheduling import mission as _mission

from ligo.skymap.io import read_sky_map
from ligo.skymap.bayestar import rasterize
from ligo.skymap import plot
from ligo.skymap.postprocess import find_greedy_credible_levels

import os, sys

from glob import glob

## command line inputs

datadir = sys.argv[1]
outdir = sys.argv[2]
allsky = sys.argv[3]
#band = sys.argv[4]
source_mag = float(sys.argv[5])

## function to estimatie the apparent AB magnitude of a counterpart
## using the 170817 U magnitude from https://arxiv.org/pdf/1710.05443.pdf for the absolute AB mag
## and getting a distance modulus from the mean GW distance estimate DISTMEAN (given in Mpc)
def estimate_apparent_ABmag(distmean,M_AB=-12.1):
    distmod = 5*np.log10(distmean*1e6) - 5
    m_AB = M_AB + distmod
    return m_AB

## function to get exposure time in a healpix pixel
def get_pixel_texp(ipix,nside,obstime,source,band):
    '''
    Inputs:
        ipix : healpix pixel index
        nside : healpix map nside
        tobs : observation time
        souce : estimated observed source rate of target
    Returns:
        pix_texp : array of [ipix,required exposure time in pixel to reach SNR of 5]
    '''
    ## get pixel sky coords
    theta, phi = hp.pix2ang(nside, ipix, nest=True)
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)
    radec = SkyCoord(ICRS(ra=ra*u.deg, dec=dec*u.deg),obstime=obstime)  
    gal = radec.galactic
    
    ## get sky background
    zodi_spec = zodi.zodi_spec_coords(radec, obstime, diag=False)
    if band=='nuv':
        galactic_spec = galactic.galactic_nuv_spec(gal.b)
        nuv_zodi = Observation(zodi_spec, nuv_band)
        nuv_galactic = Observation(galactic_spec, nuv_band)
        sky = nuv_zodi.countrate(area=area) + nuv_galactic.countrate(area=area)
    elif band=='fuv':
        galactic_spec_fuv = galactic.galactic_fuv_spec(gal.b)
        fuv_zodi = Observation(zodi_spec, fuv_band, force='extrap')
        fuv_galactic = Observation(galactic_spec_fuv, fuv_band, force='extrap')
        fuv_sky = fuv_zodi.countrate(area=area) + fuv_galactic.countrate(area=area)

    ## fuv<->nuv switch
    if band=='nuv':
        texp = uvex.sensitivity.get_exposure(source, sky)
    elif band=='fuv':
        texp = uvex.sensitivity.get_exposure(source, fuv_sky)
    else:
        print("This shouldn't happen")
        import pdb; pdb.set_trace()
    pix_texp = [ipix,texp]
    return pix_texp

## Function to get the mean exposure time across the 90% c.l. region
## with optional stats (variance, max difference from mean)
def get_max_texp(cls,nside,start_time,m_obs,band,verbose=True):
    '''
    Inputs:
        cls : 90% c.l. pixels as produced by find_greedy_credible_levels
        nside : HEALPix Nside
        start_time : observation time
        m_obs : apparent AB magnitude
        return_stats : if True, also return variance and maximum deviation from the mean exposure time
        verbose : if True, print commentary/stats
    
    Returns:
        mean_texp : mean exposure time across the 90% c.l. region
        std_texp : stadard deviation of exposure time across 90% c.l. region (optional)
        maxdiff_texp : maximum deviation from mean exposure time across 90% c.l. region (optional)   
    '''
    ## simulate observation
    sp = SourceSpectrum(ConstFlux1D, amplitude=m_obs*u.ABmag)
    obs_nuv = Observation(sp, nuv_band)
    obs_fuv = Observation(sp, fuv_band)
    source_nuv = obs_nuv.countrate(area=area)
    source_fuv = obs_fuv.countrate(area=area)
    if band=='nuv':
        source = source_nuv
    elif band=='fuv':
        source = source_fuv
    texp_list = []
    count = 0
    for ipix, cl in enumerate(cls):
        if cl <=0.9:
            pix_texp = get_pixel_texp(ipix,nside,start_time,source,band)
            if np.any(pix_texp == None):
                count += 1
                continue
            
            texp_list.append(pix_texp[1])
    if len(texp_list)==0:
        print("Warning: fits file has no valid pixels.")
        return None
    max_texp = np.max(texp_list)
    
    if verbose==True:
        print("Failed calls:", count)
        print("Number of pixels in 90%c.l. region is",len(texp_list))
        print("Maximum exposure time across 90% localization region is",max_texp)

    return max_texp


## load events
events = pd.read_csv(allsky,delimiter=',') #,skiprows=1)

## reduce to <100 sq deg localization
events_100sqdeg = events[events['area(90)'] <= 100]
ids = events_100sqdeg['simulation_id']
mags = estimate_apparent_ABmag(events_100sqdeg['distmean'],M_AB=source_mag)
events_texp = pd.DataFrame({'event_id':ids,'apparent AB mag':mags,'area(90)':events_100sqdeg['area(90)']})

## manually define some arguments
duration = u.Quantity('3 hr')
time_step = u.Quantity('1 min')
mission = getattr(_mission, 'uvex')
nside = 64
nside_hires = nside*4
delay = 0

healpix = HEALPix(nside, order='nested', frame=ICRS())
healpix_hires = HEALPix(nside_hires, order='nested', frame=ICRS())

# ## use local test events (change later to use ALL)
#local_evs = [str(i) for i in [4877,3584, 3977, 5285, 5673, 5958, 7088]]

## look at every .fits in directory
fitsfiles = glob(datadir+'/*.fits')
fitsnames = np.array(list(map(lambda filepath: filepath.split('/')[-1],fitsfiles)))

rows = []
ntot = len(ids)
progress = 1
for file_id in ids:
#for file_id in local_evs:
    ev_name = str(file_id)+'.fits'
    print("Processing",ev_name,"; (",progress,"/",ntot,")")
    progress += 1
    if ev_name not in fitsnames:
        print('Warning: event ID missing corresponding .fits for ID#',file_id)
        continue
    fitsfile = np.array(fitsfiles)[fitsnames==ev_name]
    if len(fitsfile)!=1:
        print('Warning: duplicate events with ID',file_id,'; Using first event.')
    fitsfile = fitsfile[0]    

    # Read multi-order sky map and rasterize to working resolution
    start_time = Time(fits.getval(fitsfile, 'DATE-OBS', ext=1))
    skymap_base = read_sky_map(fitsfile, moc=True)['UNIQ', 'PROBDENSITY']
    skymap = rasterize(skymap_base, healpix.level)['PROB']

    cls = find_greedy_credible_levels(skymap)

    m_obs = events_texp['apparent AB mag'][events_texp['event_id']==int(file_id)].to_numpy()[0]

    texp_output = get_max_texp(cls,nside,start_time,m_obs,band,verbose=False)
    if texp_output==None:
        print("Event",file_id,"has an invalid skymap. Increasing NSIDE...")
        skymap_hires = rasterize(skymap_base, healpix_hires.level)['PROB']
        cls_hires = find_greedy_credible_levels(skymap_hires)
        texp_output_hires = get_max_texp(cls_hires,nside_hires,start_time,m_obs,band,verbose=False)
        if texp_output_hires==None:
            print("Event",file_id,"still has an invalid skymap after increasing NSIDE. Skipping...")
            continue
        else:
            texp_max = texp_output_hires
    else:
        texp_max = texp_output
    ## create row for dataframe
    area90 = events_texp['area(90)'][events_texp['event_id']==int(file_id)].to_numpy()[0]
    rows.append([file_id,m_obs,area90,texp_max])
    print("Row added:",file_id,m_obs,area90,texp_max)
    
events_texp_stats = pd.DataFrame(rows,columns=['event_id','apparent AB mag','area(90)','texp_max (s)'])

# print(events_texp_stats)

outname = str(allsky).split('_')[-1].replace('.txt','')

events_texp_stats.to_csv(outdir+'/allsky_texp_max_'+band+'_'+outname+'.txt',index=False,sep=' ')

max_texp = 10800 ## old default 10000
events_texp_stats_10k = events_texp_stats[events_texp_stats['texp_max (s)'] <= max_texp]

# print(events_texp_stats_10k)

events_texp_stats_10k.to_csv(outdir+'/allsky_texp_max_10kscut_'+band+'_'+outname+'.txt',index=False,sep=' ')
