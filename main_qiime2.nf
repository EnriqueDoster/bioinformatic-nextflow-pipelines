#!/usr/bin/env nextflow

/*
# Example for Angus server
# Reminder to first activate qiime2 anaconda environment
/s/angus/index/common/tools/nextflow run main_qiime2.nf --reads "/s/angus/index/projs/projects_3_4/raw_sequence_data/16S_raw_combined/*_{1,2}.fq.gz" --classifier /media/AngusWorkspace/run_Jake/bioinformatic-nextflow-pipelines/gg-13-8-99-515-806-nb-classifier.qza --threads 10 -profile local_angus

# Latest qiime2 version that this script worked with was: qiime2-2020.11

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
classifier = params.classifier



Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
    .set { reads }


reads
    .map { name, forward, reverse -> [ forward.drop(forward.findLastIndexOf{"/"})[0], forward, reverse ] } //extract file name
    .map { name, forward, reverse -> [ name.toString().take(name.toString().indexOf("_")), forward, reverse ] } //extract sample name
    .map { name, forward, reverse -> [ name +"\t"+ forward +"\t"+ reverse ] } //prepare basic synthax
    .flatten()
    .collectFile(name: 'manifest.txt', newLine: true, storeDir: "${params.output}/demux", seed: "sample-id\tforward-absolute-filepath\treverse-absolute-filepath")
    .set { ch_manifest }


process Qiime2InitialProcessing {
    tag { sample_id }

    publishDir "${params.output}/read_data/", mode: "copy"


    input:
        file(manifest) from ch_manifest

    output:
        file("demux.qza") into (sample_qza)
        file("demux.qzv") into (visualization_qzv)

    """
    ${QIIME2} tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path ${manifest} \
      --output-path demux.qza \
      --input-format PairedEndFastqManifestPhred33V2

    ${QIIME2} demux summarize \
      --i-data demux.qza \
      --o-visualization demux.qzv

    """
}



process Qiime2TaxaClassification {
    tag { sample_id }

    publishDir "${params.output}/Qiime2_results", mode: "copy"

    input:
        file(qza) from sample_qza

    output:
        file("Qiime2_results/*") into (all_results)
        file("dada_results/*") into (dada_results)
    """
    ${QIIME2} dada2 denoise-paired --i-demultiplexed-seqs ${qza} \
     --o-table dada-table.qza \
     --o-representative-sequences rep-seqs.qza --p-trim-left-f ${leading) --p-trim-left-r ${trailing) --p-trunc-len-f 230 --p-trunc-len-r 230 --p-n-threads ${threads} --verbose \
     --output-dir dada_results

    ${QIIME2} phylogeny align-to-tree-mafft-fasttree \
      --i-sequences rep-seqs.qza \
      --o-alignment aligned-rep-seqs.qza \
      --o-masked-alignment masked-aligned-rep-seqs.qza \
      --o-tree unrooted-tree.qza \
      --o-rooted-tree rooted-tree.qza

    ${QIIME2} feature-classifier classify-sklearn \
      --i-classifier ${classifier} \
      --i-reads rep-seqs.qza \
      --o-classification taxonomy.qza

    ${QIIME2} taxa filter-table \
      --i-table dada-table.qza \
      --i-taxonomy taxonomy.qza \
      --p-exclude mitochondria,chloroplast \
      --o-filtered-table filtered-dada-table.qza

    unzip filtered-dada-table.qza -d exported-qiime2/
    unzip aligned-rep-seqs.qza -d exported-qiime2/
    unzip rooted-tree.qza -d exported-qiime2/
    unzip taxonomy.qza -d exported-qiime2/

    mkdir Qiime2_results/
    mv exported-qiime2/*/data/*  Qiime2_results/
    ${BIOM} convert -i Qiime2_results/feature-table.biom -o Qiime2_results/otu_table_json.biom --table-type="OTU table" --to-json

    rm demux.qza
    rm -rf exported-qiime2/
    """
}



/* Explore the addition of picrust2
qiime picrust2 full-pipeline \
   --i-table Proj_3_4_CB-dada-table-filtered.qza \
   --i-seq Proj_3_4_CB-rep-seqs.qza \
   --output-dir projs_3_4_CB-picrust2_output \
   --p-threads 1 \
   --p-hsp-method pic \
   --p-max-nsti 2 \
   --verbose
*/



/* Other visualization commands for qiime2 that could be added into the pipeline
#qiime feature-table summarize --i-table dada-table.qza \
#  --m-sample-metadata-file ${metadata} \
#  --o-visualization dada-table.qzv

#qiime feature-table tabulate-seqs --i-data rep-seqs.qza \
#  --o-visualization rep-seqs.qzv

#qiime metadata tabulate \
#  --m-input-file taxonomy.qza \
#  --o-visualization taxonomy.qzv

#qiime taxa barplot \
#  --i-table filtered-dada-table.qza \
#  --i-taxonomy taxonomy.qza \
#  --m-metadata-file ${metadata} \
#  --o-visualization filtered_taxa-bar-plots.qzv
*/




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
