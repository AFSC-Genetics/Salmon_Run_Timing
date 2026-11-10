#!/bin/bash

#SBATCH --job-name=fet_test
#SBATCH --time=2-16:00:00
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/counts_to_fet_for_localscore.%j.out
#SBATCH --partition=standard,medmem

# set variables
PROJECTDIR="/share2/nhowe/runtiming" # do not include lagging /
PREFIX="tocoho"
SUFFIX="minInd0.3_minDepthHalf_HWE"

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

Rscript --vanilla '/home/nhowe/runtiming/tocoho/scripts/localscore_FET_dumpCounts3_MM3_share2.R' ${PROJECTDIR} ${PREFIX} ${SUFFIX}

### END OF SCRIPT


