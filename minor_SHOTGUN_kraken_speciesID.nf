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

species = "Salmonella enterica"
threads = params.threads
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
Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
    .set { reads }

process RunKraken {
    tag { sample_id }
    publishDir "${params.output}/RunKraken", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".kraken.raw") > 0) "Standard/$filename"
            else if(filename.indexOf(".kraken.report") > 0) "Standard_report/$filename"
            else if(filename.indexOf(".kraken.filtered.report") > 0) "Filtered_report/$filename"
            else if(filename.indexOf(".kraken.filtered.raw") > 0) "Filtered/$filename"
            else {}
        }

    input:
       set sample_id, file(forward), file(reverse) from reads

    output:
       file("${sample_id}.kraken.report") into (kraken_report,kraken_extract_taxa)
       set sample_id, file("${sample_id}.kraken.raw") into kraken_raw
       file("${sample_id}.kraken.filtered.report") into kraken_filter_report
       file("${sample_id}.kraken.filtered.raw") into kraken_filter_raw
       file("${sample_id}.copy.R1.fastq") into forward_reads
       file("${sample_id}.copy.R2.fastq") into reverse_reads     


    """
    kraken2 --preload --db ${kraken_db} --paired ${forward} ${reverse} --threads ${threads} --report ${sample_id}.kraken.report > ${sample_id}.kraken.raw
    kraken2 --preload --db ${kraken_db} --confidence 1 --paired ${forward} ${reverse} --threads ${threads} --report ${sample_id}.kraken.filtered.report > ${sample_id}.kraken.filtered.raw
    gunzip -c ${forward} > ${sample_id}.copy.R1.fastq
    gunzip -c ${reverse} > ${sample_id}.copy.R2.fastq    
    """
}

kraken_report.toSortedList().set { kraken_l_to_w }
kraken_filter_report.toSortedList().set { kraken_filter_l_to_w }

process KrakenLongToWide {
    tag { sample_id }

    publishDir "${params.output}/KrakenLongToWide", mode: "copy"

    input:
        file(kraken_reports) from kraken_l_to_w

    output:
        file("kraken_analytic_matrix.csv") into kraken_master_matrix
        file("TaxaID.txt") into taxa_ID

    """
    grep 'Salmonella enterica' ${kraken_reports} | cut -s -f 5 | sort | uniq > TaxaID.txt
    mkdir ret
    python3 $baseDir/bin/kraken2_long_to_wide.py -i ${kraken_reports} -o ret
    mv ret/kraken_analytic_matrix.csv .
    """
}

process FilteredKrakenLongToWide {
    tag { sample_id }

    publishDir "${params.output}/Filtered_KrakenLongToWide", mode: "copy"

    input:
        file(kraken_reports) from kraken_filter_l_to_w

    output:
        file("filtered_kraken_analytic_matrix.csv") into filter_kraken_master_matrix
        file("TaxaID.txt") into filter_taxa_ID

    """
    grep 'Salmonella enterica' ${kraken_reports} | cut -s -f 5 | sort | uniq > TaxaID.txt
    mkdir ret
    python3 $baseDir/bin/kraken2_long_to_wide.py -i ${kraken_reports} -o ret
    mv ret/kraken_analytic_matrix.csv filtered_kraken_analytic_matrix.csv
    """
}

process ExtractTaxaReads {
    tag { sample_id }
    publishDir "${params.output}/TaxaReads", mode: 'copy'

    input:
       file(forward) from forward_reads
       file(taxa_ID) from taxa_ID
       set sample_id, file(raw_kraken) from kraken_raw
       file(filter_kraken) from kraken_filter_raw
       file(reverse) from reverse_reads
    output:
       file("${sample_id}.taxa.R1.fastq") into forward_taxa
       file("${sample_id}.taxa.R2.fastq") into reverse_taxa
       set sample_id , file("${sample_id}.filter.taxa.R1.fastq") into filter_forward_taxa
       file("${sample_id}.filter.taxa.R2.fastq") into filter_reverse_taxa
    """
    awk 'NR==FNR{row[\$0]=1; next} row[\$3] {print \$2}' ${taxa_ID} ${raw_kraken} > ${sample_id}.taxa.read.headers.txt
    python $baseDir/bin/extract_fastq_IDs.py ${forward} ${sample_id}.taxa.read.headers.txt ${sample_id}.taxa.R1.fastq
    python $baseDir/bin/extract_fastq_IDs.py ${reverse} ${sample_id}.taxa.read.headers.txt ${sample_id}.taxa.R2.fastq
    awk 'NR==FNR{row[\$0]=1; next} row[\$3] {print \$2}' ${taxa_ID} ${filter_kraken} > ${sample_id}.filter.taxa.read.headers.txt
    python $baseDir/bin/extract_fastq_IDs.py ${forward} ${sample_id}.filter.taxa.read.headers.txt ${sample_id}.filter.taxa.R1.fastq
    python $baseDir/bin/extract_fastq_IDs.py ${reverse} ${sample_id}.filter.taxa.read.headers.txt ${sample_id}.filter.taxa.R2.fastq

    """
}


process BlastTaxaReads {
    tag { sample_id }
    
    publishDir "${params.output}/BlastTaxaReads", mode: 'copy',  pattern: '*.out',
    saveAs: { filename ->
        if(filename.indexOf("taxa_reads_nr_filter.out") > 0) "Filtered_read_blast/$filename"
        else if(filename.indexOf("taxa_reads_nr.out") > 0) "NoFilter_read_blast/$filename"
        else {}
    }
    
    input:
       set sample_id, file(filter_forward) from filter_forward_taxa
       file(filter_reverse) from filter_reverse_taxa
       file(forward) from forward_taxa
       file(reverse) from reverse_taxa   

    output:
       set sample_id , file("${sample_id}_taxa_reads_nr_filter.out") into filter_blast_out
       file("${sample_id}_taxa_reads_nr.out") into blast_out
       
    """
    fq2fa --merge --filter $filter_forward $filter_reverse ${sample_id}.filter.interleavened.fasta
    blastn -db nt -query ${sample_id}.filter.interleavened.fasta -max_target_seqs 5 -out ${sample_id}_taxa_reads_nr_filter.out -remote -outfmt 6 

    fq2fa --merge --filter $forward $reverse ${sample_id}.interleavened.fasta
    blastn -db nt -query ${sample_id}.interleavened.fasta -max_target_seqs 5 -out ${sample_id}_taxa_reads_nr.out -remote    

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
