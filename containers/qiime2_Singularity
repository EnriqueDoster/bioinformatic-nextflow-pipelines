Bootstrap: docker
From: debian:jessie-slim

#Includes idba, trimmomatic, samtools, bwa, bedtools, freebayes, bbmap, vcftools, htslib, ncurses, kraken2, resistomeanalyzer, rarefactionanalyzer, SNPfinder

%environment

%post
    ## Jave install doesn't work, but can load java module from summit
    apt update \
    && apt install -y --no-install-recommends \
    build-essential ca-certificates sudo tcsh\
    git make automake autoconf openjdk-7-jre wget gzip unzip sed\
    zlib1g-dev curl libbz2-dev locales libncurses5-dev liblzma-dev libcurl4-openssl-dev software-properties-common apt-transport-https\
    python3-pip python3-docopt python3-pytest python-dev python3-dev\
    libcurl4-openssl-dev libssl-dev zlib1g-dev fonts-texgyre \
    gcc g++ gfortran libblas-dev liblapack-dev dos2unix\
    r-base-core r-recommended hmmer\
    && rm -rf /var/lib/apt/lists/*

    #Installing Anaconda 2 and Conda 4.5.11
    wget -c https://repo.continuum.io/archive/Anaconda2-5.3.0-Linux-x86_64.sh
    sh Anaconda2-5.3.0-Linux-x86_64.sh -bfp /usr/local

    # add bioconda channels
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

    #Download file for qiime2 installation
    wget https://data.qiime2.org/distro/core/qiime2-2019.7-py36-linux-conda.yml
    conda env create -n AmrPlusPlus_qiime2_env --file qiime2-2019.7-py36-linux-conda.yml

    . /usr/local/bin/activate AmrPlusPlus_qiime2_env

    # install bulk of bioinformatic tools using conda
    conda install biopython picrust
    
    # Run additional command for picrust
    download_picrust_files.py
    
    #ln -s /usr/local/envs/AmrPlusPlus_qiime2_env/bin/* /usr/local/bin/
    
    #Still experimenting with how to change $PATH location. 
    echo 'export PATH=$PATH:/usr/local/envs/AmrPlusPlus_qiime2_env/bin/' >> $SINGULARITY_ENVIRONMENT

%test
    #which bwa
    #which bedtools
    #which samtools
    #which freebayes
