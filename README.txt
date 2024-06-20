uvex-followup provides a set of scripts and a workflow to determine prospects for follow-up observations of gravitational-wave events with UVEX, for different target-of-opportunity selection strategies.


#################################
##      First-time setup:      ##
#################################
Clone this repository to your destination of choice:
> git clone https://github.com/criswellalexander/uvex-followup

Create a conda environment. If you decide to use a different environment name, you will need to modify the cluster submission scripts (sub_*.slurm) accordingly.
> conda create --name uvex-followup-env python=3.11

Activate the environment:
> conda activate uvex-followup-env

Install ligo.skymap via conda:
> conda install -y ligo.skymap --channel conda-forge

Install dorado-scheduling.
> pip install dorado-scheduling

Follow the instructions in the dorado-scheduling quickstart guide for CPLEX access: https://dorado-scheduling.readthedocs.io/en/latest/quickstart.html#to-set-up-the-cplex-optimization-engine

Install uvex-mission. Note that this is currently a private repo; you will need to be granted access.
The git clone command below will also ask for a password; this is NOT your Github password; you will need to create a Personal Access Token for the uvex-mission repo.
(See https://stackoverflow.com/questions/2505096/clone-a-private-repository-github, https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens)
> git clone https://github.com/uvex-mission/uvex-mission.git
Follow the installation instructions at https://github.com/uvex-mission/uvex-mission; i.e.,
> conda install -y synphot --channel conda-forge
> conda install -y specutils --channel conda-forge
> conda install -y photutils --channel conda-forge
(you should have the other dependencies already)
Then, in the uvex-mission directory:
> pip install -e .

If you intend to use the cluster submission scripts, you will also need to modify the email settings in each to
#SBATCH --mail-user=<your-email@umn.edu>
(please don't send me emails when your jobs fail)


#################################
##   uvex-followup workflow    ##
#################################

Each step can be run via a provided bash script, or you can use the various python files individually.
The workflow bash scripts are designed to work with the Minnesota Supercomputing Institute clustes' Slurm job submission framework.
You will need to adapt them slightly to your own use case.
The workflow bash scripts generally cover one [preprocess + cluster submission] segment, as you will have to wait for cluster jobs to finish running in between steps.

All parameters of interest can be modified in the params file. It is recommended to copy default_params.ini before making changes.

0 - activate your conda environment
> conda activate uvex-followup-env

1 - run workflow-part1-t_exp.sh
This will downselect to events with less than the specified localization area and submit the exposure time calculations, batched, to the cluster.
Usage: ./workflow-part1-t_exp.sh [/path/to/params_file.ini]
Example command:
> ./uvex-followup/workflow-part1-t_exp.sh ./uvex-followup/default_params.ini

2 - run workflow-part2-schedule.sh
This will ensure all exposure times are at least the specififed minimum exposure time and reformat/rebatch for scheduler submission.
Usage: ./workflow-part2-schedule.sh [/path/to/params_file.ini]
Example command:
> ./uvex-followup/workflow-part2-schedule.sh ./uvex-followup/default_params.ini

3 - run workflow-part3-postprocess
This will compute UVEX follow-up coverage for all downselected events, create plots, calculate statistics, etc..
Usage: ./workflow-part3-postprocess.sh [/path/to/params_file.ini]
Example command:
> ./uvex-followup/workflow-part3-postprocess.sh ./uvex-followup/default_params.ini
