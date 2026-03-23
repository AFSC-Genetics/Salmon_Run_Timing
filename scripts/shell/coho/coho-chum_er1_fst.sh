#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=esr1_coho-chum
#SBATCH --output=/home/nhowe/runtiming/coho/job_outfiles/coho-chum_esr1_fst_%A.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933
module load bio/angsd/0.933

chrom="NC_068428.1"

# minind 30%N, mindepth 0.5xN, max depth 10xN

angsd -b /home/nhowe/runtiming/coho/coho-chum_EE_esr1_bams.txt \
	-r ${chrom}: \
	-sites /home/nhowe/runtiming/coho/gls/coho-chum_wholegenome_minInd0.3_minDepthHalf.sites \
	-ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-out /home/nhowe/runtiming/coho/diversity/coho-chum_${chrom}_EE_minInd0.3_minDepthHalf \
	-nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -only_proper_pairs 1 \
	-C 50 -minMapQ 20 \
	-minInd 11 -setminDepth 19 -setmaxDepth 390 \
	-GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -doCounts 1

angsd -b /home/nhowe/runtiming/coho/coho-chum_LL_esr1_bams.txt \
        -r ${chrom}: \
        -sites /home/nhowe/runtiming/coho/gls/coho-chum_wholegenome_minInd0.3_minDepthHalf.sites \
        -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
        -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-out /home/nhowe/runtiming/coho/diversity/coho-chum_${chrom}_LL_minInd0.3_minDepthHalf \
        -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -only_proper_pairs 1 \
        -C 50 -minMapQ 20 \
        -minInd 11 -setminDepth 18 -setmaxDepth 370 \
        -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -doCounts 1

realSFS /home/nhowe/runtiming/coho/diversity/coho-chum_${chrom}_EE_minInd0.3_minDepthHalf.saf.idx -fold 1 /home/nhowe/runtiming/coho/diversity/coho-chum_${chrom}_LL_minInd0.3_minDepthHalf.saf.idx -fold 1 > /home/nhowe/runtiming/coho/fst/coho-chum_${chrom}_EE-LL_minInd0.3_minDepthHalf.sfs
realSFS fst index /home/nhowe/runtiming/coho/diversity/coho-chum_${chrom}_EE_minInd0.3_minDepthHalf.saf.idx -fold 1 /home/nhowe/runtiming/coho/diversity/coho-chum_${chrom}_LL_minInd0.3_minDepthHalf.saf.idx -fold 1 -sfs /home/nhowe/runtiming/coho/fst/coho-chum_${chrom}_EE-LL_minInd0.3_minDepthHalf.sfs -fstout /home/nhowe/runtiming/coho/fst/coho-chum_${chrom}_EE-LL_minInd0.3_minDepthHalf.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/coho/fst/coho-chum_${chrom}_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/coho/fst/coho-chum_${chrom}_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt

