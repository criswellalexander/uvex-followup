uvex-followup workflow

Each step can be run via a provided bash script, or you can use the various python files individually.
The workflow bash scripts generally cover one [preprocess + cluster submission] segment, as you will have to wait for cluster jobs to finish running in between steps.

1 - run workflow-part1-t_exp.sh
This will downselect to events with <100 sq. deg. localizations and submit the exposure time calculations, batched, to the cluster.
Usage: ./workflow-part1-t_exp.sh [allsky file] [skymaps dir] [out dir]
Example command:
./workflow-part1-t_exp.sh /home/vuk/crisw015/UVEX/new-sims/O5/allsky.dat /home/vuk/crisw015/UVEX/new-sims/O5/allsky/ /home/vuk/crisw015/UVEX/[OUTPUT_DIRECTORY]

2 - run workflow-part2-schedule.sh
This will ensure all exposure times are at least 500s and reformat/rebatch for scheduler submission.
Usage: ./workflow-part2-schedule.sh [skymaps dir] [out dir]
Example command:
./workflow-part2-schedule.sh /home/vuk/crisw015/UVEX/new-sims/O5/allsky/ /home/vuk/crisw015/UVEX/[OUTPUT_DIRECTORY]/

3 - run workflow-part3-postprocess
This will compute UVEX follow-up coverage for all downselected events, create plots, calculate statistics, etc..
Usage: ./workflow-part3-postprocess.sh [skymaps dir] [out dir]
Example command:
./workflow-part3-postprocess.sh /home/vuk/crisw015/UVEX/new-sims/O5/allsky/ /home/vuk/crisw015/UVEX/[OUTPUT_DIRECTORY]/


If you wish to change:

- the total tiling time considered, edit "--duration" in line 31 of sub_uvex_scheduler.slurm.
- if total tiling time is dramatically increased, consider tweaking the maximum distance in line 56 of make-coverage-plots to better reflect the increased depth of the search.
- the area cut, edit line 187 of max-texp-by-sky-loc.py and line 56 of make-coverage-plots accordingly. 
- the minimum exposure time, edit the value of min_texp on line 21 of texp_cut_and_batch.py
- the maximum exposure time, edit max_texp on line 261 of max-texp-by-sky-loc.py and line 66 of make-coverage-plots.py accordingly.
- the BNS merger rate (adjusting by some factor), edit RATE_ADJ on line 12 of workflow-part3-postprocess.sh

