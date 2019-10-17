#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

if (params.help ) {
    return help()
}
if( params.host_index ) {
    host_index = Channel.fromPath(params.host_index).toSortedList()
    //if( host_index.isEmpty() ) return index_error(host_index)
}
if( params.host ) {
    host = file(params.host)
    if( !host.exists() ) return host_error(host)
}
if( params.amr ) {
    amr = file(params.amr)
    if( !amr.exists() ) return amr_error(amr)
}
if( params.adapters ) {
    adapters = file(params.adapters)
    if( !adapters.exists() ) return adapter_error(adapters)
}
if( params.annotation ) {
    annotation = file(params.annotation)
    if( !annotation.exists() ) return annotation_error(annotation)
}
if( params.snp_annotation  ) {
    snp_annotation = file(params.snp_annotation)
}
if( params.hmm_analysis_script ) {
    hmm_analysis_script = file(params.hmm_analysis_script)
}
if( params.hmm_snp_annotation ) {
    hmm_snp_annotation = file(params.hmm_snp_annotation)
}
if( params.hmm_annotation ) {
    hmm_annotation = file(params.hmm_annotation)
}
if( params.hmm_group1 ) {
    hmm_group1 = file(params.hmm_group1)
}
if( params.hmm_group2 ) {
    hmm_group2 = file(params.hmm_group2)
}
if( params.hmm_group3 ) {
    hmm_group3 = file(params.hmm_group3)
}


threads = params.threads
smem_threads = params.smem_threads
threshold = params.threshold

min = params.min
max = params.max
skip = params.skip
samples = params.samples

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen


Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
    .set { reads }



if( !params.amr_index ) {
   process BuildAMRIndex {
       tag { amr.baseName }

       input:
           file(amr)

       output:
           file '*' into (amr_index)

       """
       bwa index ${amr}
       """
   }
}

process AlignToAMR {
     tag { sample_id }

     publishDir "${params.output}/AlignToAMR", mode: "copy"

     input:
         set sample_id, file(forward), file(reverse) from reads
         file index from amr_index
         file amr

     output:
         set sample_id, file("${sample_id}.amr.alignment.sam") into (megares_resistome_sam, megares_rarefaction_sam, megares_snp_sam , megares_snpfinder_sam)
         set sample_id, file("${sample_id}.amr.alignment.dedup.sam") into (megares_dedup_resistome_sam)
         set sample_id, file("${sample_id}.amr.alignment.dedup.bam") into (megares_dedup_resistome_bam)


     """
     ${BWA} mem ${amr} ${forward} ${reverse} -t ${threads} -R '@RG\\tID:${sample_id}\\tSM:${sample_id}' > ${sample_id}.amr.alignment.sam
     ${SAMTOOLS} view -S -b ${sample_id}.amr.alignment.sam > ${sample_id}.amr.alignment.bam
     ${SAMTOOLS} sort -n ${sample_id}.amr.alignment.bam -o ${sample_id}.amr.alignment.sorted.bam
     ${SAMTOOLS} fixmate ${sample_id}.amr.alignment.sorted.bam ${sample_id}.amr.alignment.sorted.fix.bam
     ${SAMTOOLS} sort ${sample_id}.amr.alignment.sorted.fix.bam -o ${sample_id}.amr.alignment.sorted.fix.sorted.bam
     ${SAMTOOLS} rmdup -S ${sample_id}.amr.alignment.sorted.fix.sorted.bam ${sample_id}.amr.alignment.dedup.bam
     ${SAMTOOLS} view -h -o ${sample_id}.amr.alignment.dedup.sam ${sample_id}.amr.alignment.dedup.bam
     #rm ${sample_id}.amr.alignment.bam
     #rm ${sample_id}.amr.alignment.sorted*.bam
     """
}

process RunResistome {
    tag { sample_id }

    publishDir "${params.output}/RunResistome", mode: "copy"

    input:
        set sample_id, file(sam) from megares_resistome_sam
        file annotation
        file amr

    output:
        file("${sample_id}.gene.tsv") into (megares_resistome_counts, SNP_confirm_long)

    """
    ${RESISTOME} -ref_fp ${amr} \
      -annot_fp ${annotation} \
      -sam_fp ${sam} \
      -gene_fp ${sample_id}.gene.tsv \
      -group_fp ${sample_id}.group.tsv \
      -class_fp ${sample_id}.class.tsv \
      -mech_fp ${sample_id}.mechanism.tsv \
      -t ${threshold}
    """
}

megares_resistome_counts.toSortedList().set { megares_amr_l_to_w }

process ResistomeResults {
    tag { }

    publishDir "${params.output}/ResistomeResults", mode: "copy"

    input:
        file(resistomes) from megares_amr_l_to_w

    output:
        file("AMR_analytic_matrix.csv") into amr_master_matrix

    """
    ${PYTHON3} $baseDir/bin/amr_long_to_wide.py -i ${resistomes} -o AMR_analytic_matrix.csv
    """
}




def nextflow_version_error() {
    println ""
    println "This workflow requires Nextflow version 0.25 or greater -- You are running version $nextflow.version"
    println "Run ./nextflow self-update to update Nextflow to the latest available version."
    println ""
    return 1
}

def adapter_error(def input) {
    println ""
    println "[params.adapters] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def amr_error(def input) {
    println ""
    println "[params.amr] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def annotation_error(def input) {
    println ""
    println "[params.annotation] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def fastq_error(def input) {
    println ""
    println "[params.reads] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def host_error(def input) {
    println ""
    println "[params.host] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def index_error(def input) {
    println ""
    println "[params.host_index] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def help() {
    println ""
    println "Program: AmrPlusPlus"
    println "Documentation: https://github.com/colostatemeg/amrplusplus/blob/master/README.md"
    println "Contact: Christopher Dean <cdean11@colostate.edu>"
    println ""
    println "Usage:    nextflow run main.nf [options]"
    println ""
    println "Input/output options:"
    println ""
    println "    --reads         STR      path to FASTQ formatted input sequences"
    println "    --adapters      STR      path to FASTA formatted adapter sequences"
    println "    --host          STR      path to FASTA formatted host genome"
    println "    --host_index    STR      path to BWA generated index files"
    println "    --amr           STR      path to AMR resistance database"
    println "    --annotation    STR      path to AMR annotation file"
    println "    --output        STR      directory to write process outputs to"
    println "    --KRAKENDB      STR      path to kraken database"
    println ""
    println "Trimming options:"
    println ""
    println "    --leading       INT      cut bases off the start of a read, if below a threshold quality"
    println "    --minlen        INT      drop the read if it is below a specified length"
    println "    --slidingwindow INT      perform sw trimming, cutting once the average quality within the window falls below a threshold"
    println "    --trailing      INT      cut bases off the end of a read, if below a threshold quality"
    println ""
    println "Algorithm options:"
    println ""
    println "    --threads       INT      number of threads to use for each process"
    println "    --threshold     INT      gene fraction threshold"
    println "    --min           INT      starting sample level"
    println "    --max           INT      ending sample level"
    println "    --samples       INT      number of sampling iterations to perform"
    println "    --skip          INT      number of levels to skip"
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}
