#!/bin/bash

#SBATCH --job-name=esr1PCA
#SBATCH --time=0-1:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov # update your email
#SBATCH --output=/home/nhowe/runtiming/anvil/job_outfiles/esr1_pca.%j.out # update your out file directory

module unload bio/pcangsd/0.99
module load bio/pcangsd/0.99
source /opt/bioinformatics/venv/pcangsd-0.99/bin/activate

path=/home/nhowe/runtiming/anvil
prefix=euclide
chrom=NC_068428.1

# region specific
firstpos=23141735
lastpos=23906114
gene=esr1

# call in beagle file for pop and chromosome
BEAGLE=${path}/gls/${prefix}_wholegenome_minInd0.3.beagle.gz

# cut to regions of interest
zgrep -E "${chrom}" ${BEAGLE} | awk -v s=$firstpos -v e=$lastpos -F'[\t_]' '$3 >= s && $3 <= e' | gzip > ${path}/gls/${prefix}_${chrom}_s${firstpos}_e${lastpos}_minInd0.3_${gene}.beagle.gz

pcangsd.py -threads 2 -beagle ${path}/gls/${prefix}_${chrom}_s${firstpos}_e${lastpos}_minInd0.3_${gene}.beagle.gz -o ${path}/pca/${prefix}_${chrom}_s${firstpos}_e${lastpos}_minInd0.3_${gene}

# NARROWER REGION

firstpos=23264160
lastpos=23645630

# cut to regions of interest
zgrep -E "${chrom}" ${BEAGLE} | awk -v s=$firstpos -v e=$lastpos -F'[\t_]' '$3 >= s && $3 <= e' | gzip > ${path}/gls/${prefix}_${chrom}_s${firstpos}_e${lastpos}_minInd0.3_${gene}.beagle.gz

pcangsd.py -threads 2 -beagle ${path}/gls/${prefix}_${chrom}_s${firstpos}_e${lastpos}_minInd0.3_${gene}.beagle.gz -o ${path}/pca/${prefix}_${chrom}_s${firstpos}_e${lastpos}_minInd0.3_${gene}


