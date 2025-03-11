#!/bin/bash

#SBATCH --job-name=JOBNAME  ### named job, limit 8 chars
#SBATCH --partition=nodeviper    ### partition to run on
#SBATCH --output=/scratch/tbuchlo/sim_reduced/output/chisse_chisse/peeks/%x.%j.out   ### path for output file
#SBATCH --error=/scratch/tbuchlo/sim_reduced/output/chisse_chisse/peeks/%x.%j.err    ### path for error file
#SBATCH --chdir=/scratch/tbuchlo/sim_reduced  ### some programs use a work dir
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=4 ### number of cpus assigned to job if ntasks=1
#SBATCH --mem=2GB   ### memory available for the script
#SBATCH --time=72:00:00   ### time limit (hh:mm:ss) (max is 72hr for general queue))
#SBATCH --mail-user=tbuchlo@clemson.edu    ### email to receive script notifications
#SBATCH --mail-type=BEGIN,END,FAIL    ### email notification types

/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/chisse/chisse_chisse_N0.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/chisse/chisse_chisse_N1.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/chisse/chisse_chisse_N2.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/chisse/chisse_chisse_N3.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/chisse/chisse_chisse_N4.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/chisse/chisse_chisse_N5.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/chisse/chisse_chisse_N6.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/chisse/chisse_chisse_N7.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/chisse/chisse_chisse_N8.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/chisse/chisse_chisse_N9.rev

wait