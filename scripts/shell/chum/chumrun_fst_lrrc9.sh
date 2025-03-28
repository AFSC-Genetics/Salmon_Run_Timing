#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=lrrc9
#SBATCH --output=/home/nhowe/runtiming/chum/job_outfiles/chumrun_fst_lrr9_al_%j.out
#SBATCH --error=/home/nhowe/runtiming/chum/job_outfiles/chumrun_fst_lrr9_al_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

chrom=NC_068455.1

# minInd = 10
#minDepth = 0.5*N
#Q = 20

angsd -b /home/nhowe/runtiming/chum/chumrun_early_al_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/chum/gls/chumrun_wholegenome_minDepthHalf.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_early_al_minDepthHalf -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -minInd 10 -doCounts 1 -setminDepth 13 -setmaxDepth 270.0 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/chum/chumrun_late_al_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/chum/gls/chumrun_wholegenome_minDepthHalf.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_late_al_minDepthHalf -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -minInd 10 -setminDepth 27 -setmaxDepth 550.0 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

realSFS /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_early_al_minDepthHalf.saf.idx -fold 1 /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_late_al_minDepthHalf.saf.idx -fold 1 > /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_early-late_allele_lrrc9_minDepthHalf.sfs
realSFS fst index /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_early_al_minDepthHalf.saf.idx -fold 1 /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_late_al_minDepthHalf.saf.idx -fold 1 -sfs /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_early-late_allele_lrrc9_minDepthHalf.sfs -fstout /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_early-late_allele_lrrc9_minDepthHalf.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_early-late_allele_lrrc9_minDepthHalf.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_early-late_allele_lrrc9_minDepthHalf.sfs.pbs.fst.txt

