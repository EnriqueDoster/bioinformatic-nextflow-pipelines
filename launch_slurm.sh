#!/bin/bash
#SBATCH --job-name=AMRPlusPlus
#SBATCH --partition=shas
#SBATCH --ntasks=1
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --time=23:59:00
#SBATCH --export=ALL
#SBATCH --mail-user=enriquedoster@gmail.com
#SBATCH --mail-type=ALL

module purge
module load jdk/1.8.0
module load singularity/2.5.2
module load gcc/8.2.0

./nextflow run main_AmrPlusPlus_v2_withRGI.nf -resume -profile singularity_slurm \
-w work_dir --threads 20 \
--output /scratch/summit/edoster@colostate.edu/bulktank_AMRPlusPlus_v2_results --host ../mod_bos_taurus.fna \
--reads "../bulk_tank/*_R{1,2}.fastq.gz"
