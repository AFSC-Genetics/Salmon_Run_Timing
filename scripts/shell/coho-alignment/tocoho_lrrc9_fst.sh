#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=pop_tocoho
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/tocoho_lrrc9_pop-analyses_%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

# minInd = 0.3*N, minDepth = 0.5*N, maxdepth = 10*N, minmapQ = 20

chrom=NC_034187.2

angsd -b /home/nhowe/runtiming/tocoho/tocoho_lrrc9_EE_bams.txt -r ${chrom}:25000000-26000000 -sites /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.sites -ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -anc /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -out /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_EE_25-26Mb -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -minInd 15 -setminDepth 25 -setmaxDepth 510 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/tocoho/tocoho_lrrc9_LL_bams.txt -r ${chrom}:25000000-26000000 -sites /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.sites -ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -anc /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -out /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_LL_25-26Mb -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -minInd 10 -setminDepth 16 -setmaxDepth 320 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

realSFS /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_EE_25-26Mb.saf.idx -fold 1 /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_LL_25-26Mb.saf.idx -fold 1 > /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_EE-LL_25-26Mb.sfs
realSFS fst index /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_EE_25-26Mb.saf.idx -fold 1 /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_LL_25-26Mb.saf.idx -fold 1 -sfs /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_EE-LL_25-26Mb.sfs -fstout /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_EE-LL_25-26Mb.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_EE-LL_25-26Mb.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_EE-LL_25-26Mb.sfs.pbs.fst.txt

