##NECESSARY JOB SPECIFICATIONS
#BSUB -J JobExample1         #Set the job name to "JobExample1"
#BSUB -L /bin/bash           #Uses the bash login shell to initialize the job's execution environment.
#BSUB -W 1:30                #Set the wall clock limit to 1hr and 30min
#BSUB -n 1                   #Request 1 core
#BSUB -R "span[ptile=1]"     #Request 1 core per node (This can be a maximum of 16 cores)
#BSUB -R "rusage[mem=2560]"  #Request 2560MB (2.5GB) per core (If using all cores recommend maximum is 14GB)
#BSUB -M 2560                #Set the per process enforceable memory limit to 2560MB
#BSUB -o Example1Out.%J      #Send stdout/err to "Example1Out.[jobID]"

##OPTIONAL JOB SPECIFICATIONS
#BSUB -u enriquedoster@tamu.edu       #Send all emails to email_address
#BSUB -B -N                  #Send email on job begin (-B) and end (-N)

##Your Commands After This Line

module load singularity

nextflow run main_AmrPlusPlus_v2.nf -profile singularity --output test_results
