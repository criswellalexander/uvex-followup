## this script pulls all simulated events and selects those with a sky localization < max_area
## if batch=True (default), will also create files for batch submission with N_per_batch (default 20) events/batch

import numpy as np
import pandas as pd
import os, sys, configparser
import argparse


def downselect_and_batch(allsky_file,out_dir,max_area=100, batch=True,N_batch=20):
    '''
    Function to do preprocessing for UVEX follow-up of a LVK observing run scenario.
    
    Arguments
    ------------------------
    allsky_file (str) : '/path/to/allsky_file.dat'
    out_dir (str) : '/path/to/output/directory/'
    max_area (float) : Maximum sky localization area to trigger on in sq. deg.. (default 100)
    batch (bool) : Whether to make separate batch files for scheduling submission (default True)
    N_batch (int) : Number of batch files desired (default 20)
    '''
    ## make output directory
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    ## load events
    events = pd.read_csv(allsky_file,delimiter='\t',skiprows=1)
    
    ## reduce to < max_area sq deg localization
    events_cut = events[events['area(90)'] <= max_area]
    
    percent_cut = len(events_cut)/len(events)
    
    ## some info at stdout
    print("{:0.1f}% of all events in observing scenario have less than {:0.1f} sq. deg. localization.".format(percent_cut,max_area))
    print("Total events with less than {:0.1f} sq. deg. localization: {}".format(max_area,len(events_cut)))
    
    ## save full
    savepath = os.path.join(out_dir,'allsky_cut.txt')
    events_cut.to_csv(savepath,index=False,sep=' ')
    
    if batch:
        ## make subdirectory to store batch files
        batch_dir = os.path.join(out_dir,'batches')
        os.mkdir(batch_dir)
        ## load downselected (< 100 sq. deg localization, < 10 ks max t_exp) events
        events_loaded = pd.read_csv(savepath,delimiter=' ')

        ## split into smaller chunks for batch submission
        list_of_lists = np.array_split(events_loaded,N_batch)

        batchnums = range(len(list_of_lists))
        for lst, num in zip(list_of_lists,batchnums):
            filename = 'allsky_batch'+str(num)+'.txt'
            batchpath = os.path.join(batch_dir,filename)
            lst.to_csv(batchpath,index=False,sep=',')
    
    return
    









if __name__ == '__main__':
    
    ## set up argparser
    parser = argparse.ArgumentParser(description="Prepare a LVK oberving run scenario for UVEX localization calculations.")
    parser.add_argument('params', type=str, help='/path/to/params_file.ini')
    
    args = parser.parse_args()
    
    ## set up configparser
    config = configparser.ConfigParser()
    config.read(args.params)
    
    ## get info from params file
    obs_scenario_dir = config.get("params","obs_scenario")
    out_dir          = config.get("params","save_directory")
    max_area         = float(config.get("params","max_area",fallback=100))
    N_batch          = int(config.get("params","N_batch_preproc", fallback=1))
    
    ## set additional variables accordingly
    if N_batch==1:
        batch = False
    else:
        batch = True
    
    allsky_file = obs_scenario_dir+'/allsky/allsky.dat'
    
    ## run the script
    downselect_and_batch(allsky_file,out_dir,max_area=max_area,batch=batch,N_batch=N_batch)
    
    
#     parser.add_argument('allsky_file', type=str, help='/path/to/allsky_file.dat')
#     parser.add_argument('out_dir', type=str, help='/path/to/output/directory/')
#     parser.add_argument('--max_area', type=float, help='Maximum sky localization area to trigger on in sq. deg..',default=100)
#     parser.add_argument('--nobatch', action='store_true', help='Turn off batching.')
#     parser.add_argument('--N_batch', type=int, help='Number of batch files.', default=20)
    
#     args = parser.parse_args()
    
#     batch = not args.nobatch
#     
#     downselect_and_batch(args.allsky_file,args.out_dir,max_area=args.max_area,batch=batch,N_batch=args.N_batch)
    