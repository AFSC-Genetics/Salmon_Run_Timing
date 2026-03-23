#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=pop_tocoho
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/tocoho_minDepthHalf_pop-analyses_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202


#minInd 10 minDepth Half

chrom=NC_034179.2

angsd -b /home/nhowe/runtiming/tocoho/tocoho_early_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.sites -ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -anc /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -out /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_early_minDepthHalf -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -minInd 10 -setminDepth 22 -setmaxDepth 450 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/tocoho/tocoho_late_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minDepthHalf.sites -ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -anc /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -out /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_late_minDepthHalf -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -minInd 10 -setminDepth 16 -setmaxDepth 380 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

realSFS /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_early_minDepthHalf.saf.idx -fold 1 /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_late_minDepthHalf.saf.idx -fold 1 > /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_minDepthHalf.sfs
realSFS fst index /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_early_minDepthHalf.saf.idx -fold 1 /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_late_minDepthHalf.saf.idx -fold 1 -sfs /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_minDepthHalf.sfs -fstout /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_minDepthHalf.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_minDepthHalf.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_minDepthHalf.sfs.pbs.fst.txt

