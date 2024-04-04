#!/usr/bin/env python
# coding: utf-8

## This script contains some code to:
## 1. Load in and re-aggregate the results of a batched run of max-texp-by-sky-loc.py
## 2. Ensure all exposure times are at least 500s
## 3. Rebatch and create files for use by the scheduler.

## usage: texp_cut_and_batch.py [/path/to/input/directory/with/texp/files/] [/path/to/output/directory/] [number of input batches] [number of export batches] [band]

import pandas as pd
import numpy as np
import os, sys

texp_dir = sys.argv[1]
out_dir = sys.argv[2]
N_in = int(sys.argv[3])
N_out = int(sys.argv[4])
band = sys.argv[5]

min_texp = 500.0

## load downselected (< 100 sq. deg localization, < 10 ks max t_exp) events
events = pd.read_csv(texp_dir+'allsky_texp_max_10kscut_'+band+'_batch0.txt',delimiter=' ')
for i in range(1,N_in):
    next_batch = pd.read_csv(texp_dir+'allsky_texp_max_10kscut_'+band+'_batch'+str(i)+'.txt',delimiter=' ')
    events = events.append(next_batch,ignore_index=True)

## liaise calculated t_exp to desired scheduler t_exp
## i.e., ensure minimum 500s exposures, convert s to ks, reformat for use by the scheduler
texp_sched = events['texp_max (s)'].to_numpy(copy=True)
for i, (t, ev_id) in enumerate(zip(texp_sched,events['event_id'].to_list())):
    if t <= min_texp:
        texp_sched[i] = min_texp
    else:
        continue
texp_sched = texp_sched/1000
events_texp = pd.DataFrame({'event_id':events['event_id'].tolist(),'t_exp (ks)':texp_sched})

events_texp.to_csv(os.path.join(out_dir,'allsky_sched_full.txt'),index=False,sep=' ')

## batch
batch_dir = os.path.join(out_dir,'texp_sched')
os.mkdir(batch_dir)
list_of_lists = np.array_split(events_texp,N_out)
batchnums = range(len(list_of_lists))
for lst, num in zip(list_of_lists,batchnums):
    filename = os.path.join(batch_dir,'allsky_sched_batch'+str(num)+'.txt')
    lst.to_csv(filename,index=False,sep=',')

