uvex-followup workflow

Each step can be run via a provided bash script, or you can use the various python files individually.
The workflow bash scripts are designed to work with the Minnesota Supercomputing Institute clustes' Slurm job submission framework.
You will need to adapt them slightly to your own use case.
The workflow bash scripts generally cover one [preprocess + cluster submission] segment, as you will have to wait for cluster jobs to finish running in between steps.

All parameters of interest can be modified in the params file. It is recommended to copy default_params.ini before making changes.

1 - run workflow-part1-t_exp.sh
This will downselect to events with less than the specified localization area and submit the exposure time calculations, batched, to the cluster.
Usage: ./workflow-part1-t_exp.sh [/path/to/params_file.ini]
Example command:
./uvex-followup/workflow-part1-t_exp.sh ./uvex-followup/default_params.ini

2 - run workflow-part2-schedule.sh
This will ensure all exposure times are at least the specififed minimum exposure time and reformat/rebatch for scheduler submission.
Usage: ./workflow-part2-schedule.sh [/path/to/params_file.ini]
Example command:
./uvex-followup/workflow-part2-schedule.sh ./uvex-followup/default_params.ini

3 - run workflow-part3-postprocess
This will compute UVEX follow-up coverage for all downselected events, create plots, calculate statistics, etc..
Usage: ./workflow-part3-postprocess.sh [/path/to/params_file.ini]
Example command:
./uvex-followup/workflow-part3-postprocess.sh ./uvex-followup/default_params.ini