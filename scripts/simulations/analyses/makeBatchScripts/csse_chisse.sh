#!/bin/bash

#SBATCH --job-name=c_ch_1  ### named job, limit 8 chars
#SBATCH --partition=nodeviper    ### partition to run on
#SBATCH --output=/scratch/tbuchlo/sim_reduced/output/csse_chisse/peeks/%x.%j.out   ### path for output file
#SBATCH --error=/scratch/tbuchlo/sim_reduced/output/csse_chisse/peeks/%x.%j.err    ### path for error file
#SBATCH --chdir=/scratch/tbuchlo/sim_reduced  ### some programs use a work dir
#SBATCH --cpus-per-task=62 ### number of cpus assigned to job if ntasks=1
#SBATCH --mem=3GB   ### memory available for the script
#SBATCH --time=336:00:00   ### time limit (hh:mm:ss) (max is 72hr for general queue))
#SBATCH --mail-user=tbuchlo@clemson.edu    ### email to receive script notifications
#SBATCH --mail-type=BEGIN,END,FAIL    ### email notification types

module load parallel/20220522

seq 101 200 | parallel -j 15 "
sleep {%}
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/chisse/csse_chisse_{}.rev
"