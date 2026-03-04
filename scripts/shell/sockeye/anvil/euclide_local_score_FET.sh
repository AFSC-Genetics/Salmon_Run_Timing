#!/bin/bash

#SBATCH --job-name=fet_test
#SBATCH --time=2-16:00:00
#SBATCH -c 1
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --output=/home/nhowe/runtiming/anvil/job_outfiles/counts_to_fet_for_localscore.%j.out
#SBATCH --partition=standard,medmem

# set variables
PROJECTDIR="/home/nhowe/runtiming/anvil"
PREFIX="euclide"
SUFFIX="minInd0.3_MM3"

# make sure directories exist or make them
mkdir -p ${PROJECTDIR}/localscore/

# --- Mamba Environment Activation ---------------------------------------------

# The full path to your mamba/conda installation.
# This might be on a shared file system.
CONDA_BASE=$(conda info --base)

# The path to your specific mamba environment
MAMBA_ENV_NAME="rockfish"

# Source the mamba initialization script to make 'conda' and 'mamba' commands available
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Activate the mamba environment
conda activate "${MAMBA_ENV_NAME}"

# ---- Reset the R Library Path --------------------------------------------------
# use R library set in the rockfish environment (with appropriate libraries installed)
export R_LIBS_USER="/home/nhowe/.conda/envs/rockfish/lib/R/library"

# ---- RUN R SCRIPT ---------------------------------------------------------------

Rscript --vanilla ${PROJECTDIR}'/scripts/localscore_FET_dumpCounts3_MM3.R' ${PROJECTDIR} ${PREFIX} ${SUFFIX}

### END OF SCRIPT


