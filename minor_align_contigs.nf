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
if( params.meg_mobile_annot ) {
    meg_mobile_annot = file(params.meg_mobile_annot)
    if( !meg_mobile_annot.exists() ) return meg_mobile_annot(meg_mobile_annot)
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


/*
 * Create a channel for input read files
 */
params.singleEnd = false
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { reads }


process AlignDBToContigs {
     tag { sample_id }

     publishDir "${params.output}/AlignDBsToContigs", mode: "copy"

     input:
         set sample_id, file(contigs) from reads
         file amr

     output:
         set sample_id, file("${sample_id}_blat_megares.psl") into (blat_alignments)
         file("${sample_id}_assembly_stats.txt") into (assembly_stats)

     """
     sed -i 's/ /|/g' ${contigs} 
     blat ${amr} ${contigs} ${sample_id}_blat_megares.psl
     ${PYTHON3} $baseDir/bin/assembly_stats.py ${contigs} ${sample_id} 200 1000 500
     """
}


process ParseBlatResults {
     tag { sample_id }

     publishDir "${params.output}/ParseBlatResults", mode: "copy"

     input:
         set sample_id, file(results) from blat_alignments
         file annotation
         file meg_mobile_annot

     output:
         file("${sample_id}_gene_alignments.csv") into (gene_alignments)
         file("${sample_id}_gene_alignment_stats.csv") into (stats_alignments)

     """
     ${PYTHON3} $baseDir/bin/parse_blat_psl.py ${results} ${sample_id} ${annotation} ${meg_mobile_annot}
     """
}


assembly_stats.toSortedList().set { assembly_stats_list }
gene_alignments.toSortedList().set { gene_alignments_list }
stats_alignments.toSortedList().set { stats_alignments_list }


process CombineResults {
     tag { sample_id }

     publishDir "${params.output}/Results", mode: "copy"

     input:
         file assembly_stats_list
         file gene_alignments_list
         file stats_alignments_list

     output:
         file("long_contig_assembly_stats.txt") into (long_contig_assembly_stats)
         file("long_multiple_gene_alignments.csv") into (long_multiple_gene_alignments)
         file("long_gene_alignment_stats.csv") into (long_gene_alignment_stats)

     """
     cat ${assembly_stats_list} > long_contig_assembly_stats.txt
     echo "Sample,contig,contig_start_coord,contig_end_coord,gene_header,gene_start_coord,gene_end_coord,alignment_length\n" > long_multiple_gene_alignments.csv 
     cat ${gene_alignments_list} >> long_multiple_gene_alignments.csv
     cat ${stats_alignments_list} > long_gene_alignment_stats.csv
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
def meg_mobile_annot(def input) {
    println ""
    println "[params.meg_mobile_annot] fail to open: '" + input + "' : No such file or directory"
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
