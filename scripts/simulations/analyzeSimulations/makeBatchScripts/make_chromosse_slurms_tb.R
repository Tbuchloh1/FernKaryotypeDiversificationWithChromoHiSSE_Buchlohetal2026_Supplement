setwd("/path/to/working/directory")
# settings
num_jobs_per_slurm <- 10
jobs <- 1:100
num_slurms <- length(jobs) / num_jobs_per_slurm

# read templa
template <- readLines("scripts/csse/templates/csse_csse_template.sh") # csse_csse
# template <- readLines("scripts/chisse/templates/csse_chisse_template.sh") # csse_chisse

# create job scripts
for (i in 1:num_slurms) {
    
    # copy the template
    this_template <- template
    
    # change jobname
    this_template <- gsub("JOBNAME", paste0("c_c_",i), this_template) # csse_csse
    # this_template <- gsub("JOBNAME", paste0("c_ch_",i), this_template) # csse_chisse
    
    # change job numbers
    for (j in 101:(num_jobs_per_slurm+100)) {
        
        # get old num
        oldname <- paste0("N", j-1)
        
        # get new name
        newname <- paste(j + (i - 1) * num_jobs_per_slurm)
        
        # substitute
        this_template <- gsub(oldname, newname, this_template)
        
    }
    
    # write
    this_fn <- paste0("scripts/csse/csse_csse_job_",i,".sh") # csse_csse
    # this_fn <- paste0("scripts/chisse/csse_chisse_job_",i,".sh") # csse_chisse
    writeLines(this_template, con = this_fn, sep = "\n")
       
}
