#!/bin/bash

#SBATCH --job-name=JOBNAME  ### named job, limit 8 chars
#SBATCH --partition=nodeviper    ### partition to run on
#SBATCH --output=/scratch/tbuchlo/sim_reduced/output/chisse_csse/peeks/%x.%j.out   ### path for output file
#SBATCH --error=/scratch/tbuchlo/sim_reduced/output/chisse_csse/peeks/%x.%j.err    ### path for error file
#SBATCH --chdir=/scratch/tbuchlo/sim_reduced  ### some programs use a work dir
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=2 ### number of cpus assigned to job if ntasks=1
#SBATCH --mem=1GB   ### memory available for the script
#SBATCH --time=24:00:00   ### time limit (hh:mm:ss) (max is 72hr for general queue))
#SBATCH --mail-user=tbuchlo@clemson.edu    ### email to receive script notifications
#SBATCH --mail-type=BEGIN,END,FAIL    ### email notification types

/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/csse/chisse_csse_N0.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/csse/chisse_csse_N1.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/csse/chisse_csse_N2.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/csse/chisse_csse_N3.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/csse/chisse_csse_N4.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/csse/chisse_csse_N5.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/csse/chisse_csse_N6.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/csse/chisse_csse_N7.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/csse/chisse_csse_N8.rev & 
/home/tbuchlo/software/rb_mrm_tp/rb_mrm scripts/csse/chisse_csse_N9.rev

wait