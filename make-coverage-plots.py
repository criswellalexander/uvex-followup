#!/usr/bin/env python
# coding: utf-8

## This script contains some code to make various plots for the UVEX ToO follow-up of LVK observing scenarios, plus get some basic statistics of interest.

from pathlib import Path
from astropy.table import Table, join
from astropy import units as u
from matplotlib import pyplot as plt
from matplotlib  import cm
import numpy as np
import pandas as pd
import matplotlib.image as mpimg
from matplotlib.ticker import MaxNLocator
import os, sys


def get_plots_and_stats(allsky_file,coverage_file,outdir,band,mag_AB,astro_rate,max_texp,coverage_threshold=99):
    '''
    All-in-one function for now. Should rework this to have individual fx to handle plots/stats + loop in the s2a fx from my Roman work. To-do for Mon/Tues.
    '''
    ## note that this varies from simulation to simulation
    ## find yours by doing:
    ## > sqlite3 events.sqlite 
    ## > select comment from process;
    ## This will give you the simulated rate; then compute
    ## simrate_to_astrorate = (astro_rate/sim_rate)
    ## here we use sim rate = 2.71e-6 Mpc-3 yr-1 and a fiducial astro rate = 320 Gpc-3 yr-1 = 3.2e-7 Mpc-3 yr-1
    # 2.7123599515211435e-06 1 / (Mpc3 yr)
    
    ## TODO --- INCLUDE CODE FROM ROMAN STUDY TO GO FROM SIM RATE TO ASTRO RATE, GIVEN ASTRO RATE
    
    simrate_to_astrorate = 0.11797844154885782

    ## allow for astro rate adjustment
    simrate_to_astrorate = rate_adj * simrate_to_astrorate

    ## (note that in the original Petrov+22 simulations, this value was 0.09385)

    ## set coverage threshold (demand that X% of the 90% C.I. localization region be covered. Default 99%.)
#     cov_threshold = 99

    ## load events
    events_all = pd.read_csv(allsky_file,delimiter='\t',skiprows=1)
    events_sched = pd.read_csv(coverage_file,delimiter=' ')

    ## filter to events that are well-covered and account for sun exclusion
    obs_id_list = events_sched[events_sched['percent_coverage']>=cov_threshold]['event_id'].to_list()
    obs_texp_list = events_sched[events_sched['percent_coverage']>=cov_threshold]['texp_sched (ks)'].to_list()

    events_lowcov = events_all[[(ev_id in events_sched[events_sched['percent_coverage'].to_numpy() < 5]['event_id'].to_list()) for ev_id in events_all['simulation_id'].to_list()]]

    outside_FoR_id_list = events_lowcov[(events_lowcov['area(90)'].to_numpy() < 100) & (events_lowcov['distmean'].to_numpy() < 300)]['simulation_id'].to_list()


    events_texp = pd.read_csv(outdir+'/texp_out/'+'allsky_texp_max_'+band+'_batch0.txt',delimiter=' ')
    for i in range(1,40):
        next_batch = pd.read_csv(outdir+'/texp_out/'+'allsky_texp_max_'+band+'_batch'+str(i)+'.txt',delimiter=' ')
        events_texp = events_texp.append(next_batch,ignore_index=True)
    texp_cut_id_list = events_texp[events_texp['texp_max (s)'] > max_texp]['event_id'].to_list()
    events_texp_reject = events_all[[(sid in texp_cut_id_list) for sid in events_all['simulation_id']]]
    events_texp_reject_id_list = events_texp_reject['simulation_id'].to_list()

    ## anything, regardless of other factors, that gets <1% coverage gets marked with black x
    ## mark out # of tiles = 0 FoR
    ## texp long and rejected gets marked with grey x
    ## three rejection categories: completely excluded, partially excluded, can't be tiled under time constraints
    events_notcov = events_all[[(ev_id in events_sched[events_sched['percent_coverage'].to_numpy() < cov_threshold]['event_id'].to_list()) for ev_id in events_all['simulation_id'].to_list()]]
    events_0cov = events_notcov[[(ev_id in events_sched[events_sched['percent_coverage'].to_numpy() < 0.01]['event_id'].to_list()) for ev_id in events_notcov['simulation_id'].to_list()]]
    events_0cov_id_list = events_0cov['simulation_id'].to_list()
    events_semicov = events_notcov[[(ev_id in events_sched[events_sched['percent_coverage'].to_numpy() >= 0.01]['event_id'].to_list()) for ev_id in events_notcov['simulation_id'].to_list()]]
    events_semicov_sched = events_sched[[(ev_id in events_semicov['simulation_id'].to_list()) for ev_id in events_sched['event_id'].to_list()]]


    ## see if event is in principle coverable under ideal circumstances
    coverable_filt = np.floor((max_texp/1e3)/events_semicov_sched['texp_sched (ks)'].to_numpy()) <= np.ceil(events_semicov['area(90)'].to_numpy()/10) #True if coverable in perfect world
    events_coverable = events_semicov[coverable_filt]
    events_coverable_id_list = events_coverable['simulation_id'].to_list()
    events_timex = events_semicov[~coverable_filt]
    events_timex_id_list = events_timex['simulation_id'].to_list()

    ## get which events in the full simulation are observed
    ## True if observed, False if not
    obs_filter = [(event_id in obs_id_list) for event_id in events_all['simulation_id'].to_list()]
    # FoR_filter = [(event_id in outside_FoR_id_list) for event_id in events_all['simulation_id'].to_list()]
    full_exclusion_filt = [(event_id in events_0cov_id_list) for event_id in events_all['simulation_id'].to_list()]
    partial_exclusion_filt = [(event_id in events_coverable_id_list) for event_id in events_all['simulation_id'].to_list()]
    time_exclusion_filt = [(event_id in events_timex_id_list) for event_id in events_all['simulation_id'].to_list()]
    texp_reject_filt = [(event_id in events_texp_reject_id_list) for event_id in events_all['simulation_id'].to_list()]
    combined_filter = np.array(obs_filter) | np.array(full_exclusion_filt) | np.array(partial_exclusion_filt) | np.array(time_exclusion_filt)

    ## get statistics on number of events *in princple* coverable vs. those we actually observe
    ## THESE SHOULD BE UPDATED LATER TO ACCOUNT FOR THE CATEGORIZATION ABOVE
    events_midcov = events_sched[(events_sched['percent_coverage'] < cov_threshold)].copy()
    midcov_areas = [events_all[events_all['simulation_id']==ev_id]['area(90)'].to_numpy()[0] for ev_id in events_midcov['event_id']]
    events_midcov['area(90)'] = midcov_areas
    midcov_frac_coverable = len(events_midcov[np.round((max_texp/1e3)/events_midcov['texp_sched (ks)']) > np.round(events_midcov['area(90)']/10)])/len(events_midcov)
    midcov_ncoverable = midcov_frac_coverable*len(events_midcov)
    frac_obs_v_coverable = np.sum(obs_filter)/(np.sum(obs_filter) + midcov_ncoverable)

    ## generate event lists for plot
    obs_dist = events_all[obs_filter]['distmean'].to_list()
    obs_area = events_all[obs_filter]['area(90)'].to_list()
    unobs_dist = events_all[np.invert(combined_filter)]['distmean'].to_list()
    unobs_area = events_all[np.invert(combined_filter)]['area(90)'].to_list()
    full_exclusion_dist = events_all[full_exclusion_filt]['distmean'].to_list()
    full_exclusion_area = events_all[full_exclusion_filt]['area(90)'].to_list()
    partial_exclusion_dist = events_all[partial_exclusion_filt]['distmean'].to_list()
    partial_exclusion_area = events_all[partial_exclusion_filt]['area(90)'].to_list()
    time_exclusion_dist = events_all[time_exclusion_filt]['distmean'].to_list()
    time_exclusion_area = events_all[time_exclusion_filt]['area(90)'].to_list()
    texp_reject_dist = events_all[texp_reject_filt]['distmean'].to_list()
    texp_reject_area = events_all[texp_reject_filt]['area(90)'].to_list()



    # outside_FoR_dist = events_all[FoR_filter]['distmean'].to_list()
    # outside_FoR_area = events_all[FoR_filter]['area(90)'].to_list()

    ## save all info for selected events
    selected_filter = [(event_id in obs_id_list) for event_id in events_sched['event_id'].to_list()]
    events_selected_gwinfo = events_all[obs_filter].copy()
    events_selected_eminfo = events_sched[selected_filter].copy()
    events_selected_allinfo = pd.concat([events_selected_gwinfo.reset_index(drop=True),events_selected_eminfo.reset_index(drop=True)],axis=1)
    csv_savename = os.path.join(outdir,"allsky_selected_em.txt".format(band,mag_AB))
    events_selected_eminfo.to_csv(csv_savename,index=False)
    csv_savename_gwem = os.path.join(outdir,"allsky_selected_gwem.txt".format(band,mag_AB))
    events_selected_allinfo.to_csv(csv_savename_gwem,index=False)

    stat_savename = os.path.join(outdir,"uvex_event_statistics_{}_MAB{:0.1f}.txt".format(band,mag_AB))
    print("Calculating statistics...")
    with open(stat_savename,'w') as outfile:
        print("Number of selected events is",len(obs_dist),"; this is",100*len(obs_dist)/len(events_all),"percent of the catalogue.",file=outfile)
        ## covert to predictions for yearly rate
        print("Warning: using a conversion factor from simulated rate to astrophysical rate of {:0.4f}.".format(simrate_to_astrorate))
        print("Check that this is accurate for your simulations; see comments at beginning of this script for more details.")
        print("The following predictions use a (astrophysical rate / simulated rate) conversion factor of {:0.4f}.".format(simrate_to_astrorate),file=outfile)
        print("Total number of events is ",len(events_all)*simrate_to_astrorate,file=outfile)
        print("Predicted number of selected events (in 1 yr) is ",len(obs_dist)*simrate_to_astrorate,file=outfile)
        frac_obs_v_coverable = np.sum(obs_filter)/(np.sum(obs_filter) + midcov_ncoverable)
        print("Fraction of coverable events within UVEX field of regard:",frac_obs_v_coverable,file=outfile)
        print("(i.e., fraction of events lost to some level of sun exclusion is {:0.2f})".format((1-frac_obs_v_coverable)),file=outfile)
        ## exposure time statistics
        obs_texp_arr = np.array(obs_texp_list)*1000
        median_texp_sel = np.median(obs_texp_arr)
        min_texp_sel = np.min(obs_texp_arr)
        max_texp_sel = np.max(obs_texp_arr)
        print("Exposure time statistics for selected events (in s):",file=outfile)
        print("Median:", median_texp_sel,"; Min:",min_texp_sel,"; Max:",max_texp_sel,file=outfile)
        ## tiling statistics
        obs_tile_list = events_sched[[(ev_id in obs_id_list) for ev_id in events_sched['event_id']]]['tiles_to_99pct'].to_list()
        obs_tile_arr = np.array(obs_tile_list)
        median_tile_sel = np.median(obs_tile_arr)
        min_tile_sel = np.min(obs_tile_arr)
        max_tile_sel = np.max(obs_tile_arr)
        print("Tiling statistics for selected events (tiles to 99% coverage of 90% localization):",file=outfile)
        print("Median:", median_tile_sel,"; Min:",min_tile_sel,"; Max:",max_tile_sel,file=outfile)
        frac = np.sum(obs_tile_arr <=5)/len(obs_tile_arr)
        print("Fraction of selected events covered in <5 tiles is ", frac, file=outfile)
        print("This corresponds to {} predicted events in 1 yr.".format(frac*len(obs_dist)*simrate_to_astrorate),file=outfile)
    print("Statistics saved to {}.".format(stat_savename))

    print("Making plots...")
    with plt.style.context('seaborn-talk'):
        ax = plt.axes()
        ax.set_xscale('log')
        ax.set_yscale('log')
    #     ax.set_yscale('log')
        ax.set_xlim(1e2, 1e3)
        ax.set_ylim(1, 1e4)
        ax.grid(which='both', axis='x')
        ax.grid()
        # ax.plot(O5_events['distmean'], O5_events['area(90)'], '.', ms=6, label='O5 events')
        ## unselected events
        plt.scatter(unobs_dist, unobs_area, marker='.', s=8, label='Unselected events',color='darkgray')
        ## selected but unobserved by category
        plt.scatter(full_exclusion_dist+partial_exclusion_dist, full_exclusion_area+partial_exclusion_area, marker='x', s=14,linewidth=0.75, 
                label='Complete or partial field of regard exclusion', color='black')
    #     plt.scatter(partial_exclusion_dist, partial_exclusion_area, marker='x', s=10, 
    #             label='Partial field of regard exclusion', color='slategrey')
        plt.scatter(time_exclusion_dist, time_exclusion_area, marker='x', s=15,linewidth=0.75,  
                label='Total epoch time constraint', color='darkorange')
        plt.scatter(texp_reject_dist, texp_reject_area, marker='x', s=15,linewidth=0.75,  
                label='__nolabel__', color='darkorange')
        ##observed
        plt.scatter(obs_dist, obs_area, marker='o', s=17,linewidth=1, 
                    label='Selected events',c=np.array(obs_texp_list)*1000,cmap=cm.cool_r)

        #ax.fill_between(DL, A, ax1.get_ylim()[1], color='lightgray')
        #ax.plot(DL, A, color='gray', lw=4)

        ax.set_xlabel('Luminosity distance (Mpc)')
        ax.set_ylabel('90% credible area (deg$^2$)')
        cbar = plt.colorbar()
        cbar.set_label('Exposure Time (s)')#, rotation=270)
        plt.legend(loc="upper left")
        plt.title('Event Selection for UVEX ToO with {}'.format(band.upper())+' $M_{AB}=$'+'${:0.1f}$'.format(mag_AB))
        plot_savename = os.path.join(outdir,'uvex_event_selection_{}_MAB{:0.1f}.png'.format(band,mag_AB))
        plt.savefig(plot_savename,bbox_inches='tight',dpi=300)
        plt.close()

    plt.figure()
    plt.hist(obs_texp_arr,bins=20,color='mediumorchid')
    plt.title("Histogram of Exposure Times")
    plt.xlabel("Exposure Time (s)")
    plt.ylabel("Count")
    hist1_savename = os.path.join(outdir,'uvex_texp_histogram_{}_MAB{:0.1f}.png'.format(band,mag_AB))
    plt.savefig(hist1_savename,bbox_inches='tight')
    plt.close()

    plt.figure()
    ax = plt.figure().gca()
    plt.hist(obs_tile_arr,bins=20,color='mediumorchid')
    plt.title("Histogram of Tiles to Reach {}% Coverage of 90% c.l. Localization Region".format(cov_threshold))
    plt.xlabel("Tiles")
    plt.ylabel("Count")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    hist2_savename = os.path.join(outdir,'uvex_tiling_histogram_{}_MAB{:0.1f}.png'.format(band,mag_AB))
    plt.savefig(hist2_savename,bbox_inches='tight')
    plt.close()
    
    return

# allsky_file = sys.argv[1]
# coverage_file = sys.argv[2]
# outdir = sys.argv[3]
# band = sys.argv[4]
# mag_AB = float(sys.argv[5])

# if len(sys.argv)==7:
#     rate_adj = float(sys.argv[6])
# else:
#     rate_adj = 1.0
    
if __name__ == '__main__':
    
    ## set up argparser
    parser = argparse.ArgumentParser(description="Given observing schedules, compute tiling and statistics.")
    parser.add_argument('params', type=str, help='/path/to/params_file.ini')
    
    args = parser.parse_args()
    
    ## set up configparser
    config = configparser.ConfigParser()
    config.read(args.params)
    
    ## get info from params file
    obs_scenario_dir   = config.get("params","obs_scenario")
    out_dir            = config.get("params","save_directory")
    band               = config.get("params","band")
    source_mag         = float(config.get("params","KNe_mag_AB"))
    astro_rate_median  = float(config.get("params","astro_bns_median"))
    astro_rate_bounds  = eval(str(config.get("params","astro_bns_interval_90")))
    max_texp           = config.get("params","max_texp")
    coverage_threshold = config.get("params","coverage_threshold",fallback=99)
    
    allsky_file = obs_scenario_dir+'/allsky/allsky.dat'
    coverage_file = outdir+'/allsky_coverage.txt'
    
    astro_rate = [astro_rate_bounds[0], astro_rate_median, astro_rate_bounds[1]]
    
    ## run the script
    get_plots_and_stats(allsky_file,coverage_file,outdir,band,source_mag,astro_rate,max_texp,coverage_threshold)