#!/bin/bash -l
#PBS -l walltime=95:00:00,mem=5gb,nodes=1:ppn=2
#PBS -o /scratch.global/test_WGS/WGS_SNP_pipelines/outfile_o_file
#PBS -e /scratch.global/test_WGS/WGS_SNP_pipelines/errorfile_e_file
#PBS -q mesabi

#PBS -m abe
#PBS -M edoster@umn.edu

module purge
module load singularity

cd /scratch.global/test_WGS/WGS_SNP_pipelines

/home/noyes046/shared/tools/nextflow run main_combined_pipeline.nf --reference_genome /scratch.global/test_WGS/ref_L_monocytogenes_NC_003210.fasta --reads '/panfs/roc/risdb_new/genometrakr/listeria_monocytogenes/*_{1,2}.fastq.gz' -profile singularity_pbs --output /scratch.global/test_WGS/WGS_SNP_pipelines/GenomeTrakr_Listeria_qsub --threads 128 -w qsub_work -resume
