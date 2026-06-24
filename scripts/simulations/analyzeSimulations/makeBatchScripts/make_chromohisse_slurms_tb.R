setwd("/path/to/working/directory")
# settings
num_jobs_per_slurm <- 10
jobs <- 1:100
num_slurms <- length(jobs) / num_jobs_per_slurm

# read template
# template <- readLines("scripts/chisse/templates/chisse_chisse_template.sh") # chisse_chisse
template <- readLines("scripts/csse/templates/chisse_csse_template.sh") # chisse_csse

# create job scripts
for (i in 1:num_slurms) {
    
    # copy the template
    this_template <- template
    
    # change jobname
    # this_template <- gsub("JOBNAME", paste0("ch_ch_",i), this_template) # chisse_chisse
    this_template <- gsub("JOBNAME", paste0("ch_c_",i), this_template) # chisse_csse
    
    # change job numbers
    for (j in 1:num_jobs_per_slurm) {
        
        # get old num
        oldname <- paste0("N", j-1)
        
        # get new name
        newname <- paste(j + (i - 1) * num_jobs_per_slurm)
        
        # substitute
        this_template <- gsub(oldname, newname, this_template)
        
    }
    
    # write
    # this_fn <- paste0("scripts/chisse/chisse_chisse_job_",i,".sh") # chisse_chisse
    this_fn <- paste0("scripts/csse/chisse_csse_job_",i,".sh") # chisse_csse
    writeLines(this_template, con = this_fn, sep = "\n")
       
}
