#!/bin/bash

#SBATCH --job-name=esr1PCA
#SBATCH --time=0-1:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov # update your email
#SBATCH --output=/home/nhowe/runtiming/coho/job_outfiles/esr1_pca.%j.out # update your out file directory

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

# Species/run
path=/home/nhowe/runtiming/coho
prefix=coho-chum
suffix=minInd0.3_minDepthHalf

# gene region
gene=esr1
chrom=NC_068428.1

# BROADER REGION
firstpos=23141735
lastpos=23906114

# call in beagle file for pop and chromosome
BEAGLE=${path}/gls/${prefix}_wholegenome_${suffix}.beagle.gz

# cut to regions of interest
zgrep -E "${chrom}" ${BEAGLE} | awk -v s=$firstpos -v e=$lastpos -F'[\t_]' '$3 >= s && $3 <= e' | gzip > ${path}/gls/${prefix}_${chrom}_s${firstpos}_e${lastpos}_${suffix}_${gene}.beagle.gz
pcangsd.py -threads 2 -beagle ${path}/gls/${prefix}_${chrom}_s${firstpos}_e${lastpos}_${suffix}_${gene}.beagle.gz -o ${path}/pca/${prefix}_${chrom}_s${firstpos}_e${lastpos}_${suffix}_${gene}

# NARROWER REGION

firstpos=23264160
lastpos=23645630

# cut to regions of interest
zgrep -E "${chrom}" ${BEAGLE} | awk -v s=$firstpos -v e=$lastpos -F'[\t_]' '$3 >= s && $3 <= e' | gzip > ${path}/gls/${prefix}_${chrom}_s${firstpos}_e${lastpos}_${suffix}_${gene}.beagle.gz
pcangsd.py -threads 2 -beagle ${path}/gls/${prefix}_${chrom}_s${firstpos}_e${lastpos}_${suffix}_${gene}.beagle.gz -o ${path}/pca/${prefix}_${chrom}_s${firstpos}_e${lastpos}_${suffix}_${gene}
