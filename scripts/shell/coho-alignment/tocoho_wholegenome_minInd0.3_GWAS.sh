#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=gwas_tocoho
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/tocoho_wholegenome_gwas_%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=himem
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933
module load bio/angsd/0.933

# Angsd code does not work with underscore after "NC_"
# error says that it cannot find chr:NC in fai index, so it must tab on underscores. After removing underscore, it seems to work just fine
#zcat /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.beagle.gz | sed 's!NC_03!NC03!g' | gzip > /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_rmNCunderscore.beagle.gz

angsd -doMaf 4 -beagle /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_rmNCunderscore.beagle.gz -yBin /home/nhowe/runtiming/tocoho/scripts/tocoho_runtime_GWAS_bin.txt -doAsso 6 -out /home/nhowe/runtiming/tocoho/gwas/tocoho_gwas_wholegenome -fai /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic_rmNCunderscore.fna.fai
