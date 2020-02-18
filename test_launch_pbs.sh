#!/bin/bash -l
#PBS -l walltime=95:00:00,mem=5gb,nodes=1:ppn=2
#PBS -o /scratch.global/run_proj7/outfile_o_file
#PBS -e /scratch.global/run_proj7/errorfile_e_file
#PBS -q mesabi

#PBS -m abe
#PBS -M edoster@umn.edu

module purge
module load singularity

cd /scratch.global/run_proj7/amrplusplus_v2/

nextflow run minor_to_nonhost.nf -profile MSI_pbs -w /scratch.global/run_proj7/work2_dir_nonhost --threads 24 --reads '/home/noyes046/shared/projects/proj7_raw_reads/*R{1,2}.fastq.gz' --output /scratch.global/run_proj7/proj7_nonhost_fastq --host /panfs/roc/risdb/genomes/Bos_taurus/Bos_taurus_UMD_3.1/bwa/Bos_taurus_UMD_3.1.fa --adapters /panfs/roc/msisoft/trimmomatic/0.33/adapters/all_illumina_adapters.fa --kraken_db /home/noyes046/shared/databases/kraken2_databases/Rumen_kraken_v2_Nov2019/ -resume -with-report non_host_run.report -with-trace -with-timeline

