#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=pop_pink-chum
#SBATCH --output=/home/nhowe/runtiming/pink/job_outfiles/pink-chum_lrrc9_allele_fst_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

chrom=NC_068455.1

angsd -b /home/nhowe/runtiming/pink/pink-chum_EE_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/pink/gls/pink-chum_wholegenome_minInd0.3.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/pink/diversity/pink-chum_${chrom}_EE_minInd0.3 -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 15 -doCounts 1 -minInd 19 -setminDepth 65 -setmaxDepth 1300 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/pink/pink-chum_LL_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/pink/gls/pink-chum_wholegenome_minInd0.3.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/pink/diversity/pink-chum_${chrom}_LL_minInd0.3 -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 15 -doCounts 1 -minInd 14 -setminDepth 47 -setmaxDepth 940 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

realSFS /home/nhowe/runtiming/pink/diversity/pink-chum_${chrom}_EE_minInd0.3.saf.idx -fold 1 /home/nhowe/runtiming/pink/diversity/pink-chum_${chrom}_LL_minInd0.3.saf.idx -fold 1 > /home/nhowe/runtiming/pink/fst/pink-chum_${chrom}_EE-LL_minInd0.3.sfs
realSFS fst index /home/nhowe/runtiming/pink/diversity/pink-chum_${chrom}_EE_minInd0.3.saf.idx -fold 1 /home/nhowe/runtiming/pink/diversity/pink-chum_${chrom}_LL_minInd0.3.saf.idx -fold 1 -sfs /home/nhowe/runtiming/pink/fst/pink-chum_${chrom}_EE-LL_minInd0.3.sfs -fstout /home/nhowe/runtiming/pink/fst/pink-chum_${chrom}_EE-LL_minInd0.3.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/pink/fst/pink-chum_${chrom}_EE-LL_minInd0.3.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/pink/fst/pink-chum_${chrom}_EE-LL_minInd0.3.sfs.pbs.fst.txt

