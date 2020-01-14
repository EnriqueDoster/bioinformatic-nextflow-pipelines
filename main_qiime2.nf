#!/usr/bin/env nextflow

/*
qiime dada2 denoise-paired --i-demultiplexed-seqs /s/angus/index/projs/projects_3_4/analysis/Qiime2/Proj_3_4_CB/Proj_3_4_CB-paired-end.qza --o-table /s/angus/index/projs/projects_3_4/analysis/Qiime2/Proj_3_4_CB/Proj_3_4_CB-dada-table.qza --o-representative-sequences /s/angus/index/projs/projects_3_4/analysis/Qiime2/Proj_3_4_CB/Proj_3_4_CB-rep-seqs.qza --p-trim-left-f 5 --p-trim-left-r 5 --p-trunc-len-f 240 --p-trunc-len-r 240 --p-n-threads 8 --verbose

/s/angus/index/common/tools/nextflow run main_qiime2.nf --reads "/s/angus/index/projs/projects_3_4/raw_sequence_data/16S_raw_combined/*_{1,2}.fq.gz" --metadata "/media/AngusWorkspace/proj3_16S/Proj_3_metadata_full.tsv" -profile local_angus

*/


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

if(params.kraken_db) {
    kraken_db = file(params.kraken_db)
}

threads = params.threads

threshold = params.threshold

min = params.min
max = params.max
skip = params.skip
samples = params.samples

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen

metadata = params.metadata


/*
qiime picrust2 full-pipeline \
   --i-table Proj_3_4_CB-dada-table-filtered.qza \
   --i-seq Proj_3_4_CB-rep-seqs.qza \
   --output-dir projs_3_4_CB-picrust2_output \
   --p-threads 1 \
   --p-hsp-method pic \
   --p-max-nsti 2 \
   --verbose


   qiime feature-table tabulate-seqs --i-data Proj_3_4_CB-rep-seqs.qza --o-visualization Proj_3_4_CB-rep-seqs.qzv
   qiime feature-classifier classify-sklearn --i-classifier gg-13-8-99-nb-classifier.qza --i-reads Proj_3_4_CB-rep-seqs.qza --o-classification Proj_3_4_CB-taxonomy.qza
   qiime taxa filter-table --i-table Proj_3_4_CB-dada-table.qza --i-taxonomy Proj_3_4_CB-taxonomy.qza --p-exclude mitochondria,chloroplast --o-filtered-table Proj_3_4_CB-dada-table-filtered.qza
   qiime metadata tabulate --m-input-file Proj_3_4_CB-taxonomy.qza --o-visualization Proj_3_4_CB-taxonomy.qzv
   qiime feature-table summarize --i-table Proj_3_4_CB-dada-table-filtered.qza --m-sample-metadata-file Proj_3_4_CB_metadata_full.tsv --o-visualization Proj_3_4_CB-dada-table-filtered.qzv
   qiime taxa barplot --i-table Proj_3_4_CB-dada-table-filtered.qza --i-taxonomy Proj_3_4_CB-taxonomy.qza --m-metadata-file Proj_3_4_CB_metadata_full.tsv --o-visualization Proj_3_4_CB-taxonomy-filtered.qzv
   qiime alignment mafft --i-sequences Proj_3_4_CB-rep-seqs.qza --o-alignment Proj_3_4_CB-aligned.qza
   qiime alignment mask --i-alignment Proj_3_4_CB-aligned.qza --o-masked-alignment Proj_3_4_CB-aligned-masked.qza
   qiime phylogeny fasttree --i-alignment Proj_3_4_CB-aligned-masked.qza --o-tree Proj_3_4_CB-aligned-masked-unrooted.qza
   qiime phylogeny midpoint-root --i-tree Proj_3_4_CB-aligned-masked-unrooted.qza --o-rooted-tree Proj_3_4_CB-aligned-masked-rooted.qza


*/


Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
    .set { reads }


reads
    .map { name, forward, reverse -> [ forward.drop(forward.findLastIndexOf{"/"})[0], forward, reverse ] } //extract file name
    .map { name, forward, reverse -> [ name.toString().take(name.toString().indexOf("_")), forward, reverse ] } //extract sample name
    .map { name, forward, reverse -> [ name +","+ forward + ",forward\n" + name +","+ reverse +",reverse" ] } //prepare basic synthax
    .flatten()
    .collectFile(name: 'manifest.txt', newLine: true, storeDir: "${params.output}/demux", seed: "sample-id,absolute-filepath,direction")
    .set { ch_manifest }

/* This section mostly works, but the header for the manifest still contains the "\t" characters
reads
    .map { name, forward, reverse -> [ forward.drop(forward.findLastIndexOf{"/"})[0], forward, reverse ] } //extract file name
    .map { name, forward, reverse -> [ name.toString().take(name.toString().indexOf("_")), forward, reverse ] } //extract sample name
    .map { name, forward, reverse -> [ name +"\t"+ forward + "\t"+ reverse] } //prepare basic synthax
    .flatten()
    .collectFile(name: 'manifest.txt', newLine: true, storeDir: "${params.output}/demux", seed: "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
    .set { ch_manifest }
*/


process Qiime2InitialProcessing {
    tag { sample_id }

    publishDir "${params.output}/Qiime2_initial_processing", mode: "copy"


    input:
        file(manifest) file(reverse) from ch_manifest

    output:
        file("demux.qza") into (ch_qiime2_raw)

    """
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path ${manifest} \
      --output-path demux.qza \
      --input-format PairedEndFastqManifestPhred33

    qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza \
     --o-table dada-table.qza \
     --o-representative-sequences rep-seqs.qza --p-trim-left-f 5 --p-trim-left-r 5 --p-trunc-len-f 240 --p-trunc-len-r 240 --p-n-threads ${threads} --verbose

    qiime feature-table summarize --i-table dada-table.qza \
      --m-sample-metadata-file ${metadata} \
      --o-visualization dada-table.qzv

    qiime feature-table tabulate-seqs --i-data rep-seqs.qza \
      --o-visualization rep-seqs.qzv

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
