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

kraken_db = params.kraken_db

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

process RunKraken {
    tag { sample_id }
    publishDir "${params.output}/RunKraken", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".kraken.raw") > 0) "Standard/$filename"
            else if(filename.indexOf(".kraken.report") > 0) "Standard/$filename"
            else if(filename.indexOf(".kraken.filtered.report") > 0) "Filtered/$filename"
            else if(filename.indexOf(".kraken.filtered.raw") > 0) "Filtered/$filename"
            else {}
        }

    input:
       set sample_id, file(forward), file(reverse) from reads

    output:
       file("${sample_id}.kraken.filtered.report") into kraken_filter_report
       file("${sample_id}.kraken.report") into kraken_report

    """
    kraken2 --preload --db ${kraken_db} --paired ${forward} ${reverse} --threads ${threads} --report ${sample_id}.kraken.report > ${sample_id}.kraken.raw
    kraken2 --preload --db ${kraken_db} --confidence 1 --paired ${forward} ${reverse} --threads ${threads} --report ${sample_id}.kraken.filtered.report > ${sample_id}.kraken.filtered.raw
    """
}

kraken_report.toSortedList().set { kraken_l_to_w }
kraken_filter_report.toSortedList().set { kraken_filter_l_to_w }

process KrakenLongToWide {
    tag { }

    publishDir "${params.output}/KrakenLongToWide", mode: "copy"

    input:
        file(kraken_reports) from kraken_l_to_w

    output:
        file("kraken_analytic_matrix.csv") into kraken_master_matrix

    """
    mkdir ret
    python3 $baseDir/bin/kraken2_long_to_wide.py -i ${kraken_reports} -o ret
    mv ret/kraken_analytic_matrix.csv .
    """
}

process FilteredKrakenLongToWide {
    tag { }

    publishDir "${params.output}/Filtered_KrakenLongToWide", mode: "copy"

    input:
        file(kraken_reports) from kraken_l_to_w

    output:
        file("filtered_kraken_analytic_matrix.csv") into kraken_master_matrix

    """
    mkdir ret
    python3 $baseDir/bin/kraken2_long_to_wide.py -i ${kraken_reports} -o ret
    mv ret/filtered_kraken_analytic_matrix.csv .
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
