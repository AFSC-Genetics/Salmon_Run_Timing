#!/bin/bash

#SBATCH --job-name=fet_test
#SBATCH --time=2-16:00:00
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/chum/job_outfiles/counts_to_fet_for_localscore.%j.out
#SBATCH --partition=standard,medmem

# set variables
PROJECTDIR="/home/nhowe/runtiming/chum"
PREFIX="chumrun"
SUFFIX="minInd0.3_minDepthHalf_MM3"

# make sure directories exist or make them
mkdir -p ${PROJECTDIR}"/localscore/"

# --- Mamba Environment Activation ---------------------------------------------

# The full path to the mamba installation.
CONDA_BASE=$(conda info --base)

# Source the mamba initialization script to make 'conda' and 'mamba' commands available
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Activate the mamba environment
# Rockfish environment has the appropriate libraries installed
conda activate "rockfish"

# reset R library to the rockfish environment
export R_LIBS_USER="/home/nhowe/.conda/envs/rockfish/lib/R/library"

# ---- RUN R SCRIPT ---------------------------------------------------------------

Rscript --vanilla ${PROJECTDIR}'/scripts/localscore_FET_dumpCounts3_MM3.R' ${PROJECTDIR} ${PREFIX} ${SUFFIX}

### END OF SCRIPT


