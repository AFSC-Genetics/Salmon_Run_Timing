#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=sock-all_lrrc9_fst
#SBATCH --output=/home/nhowe/runtiming/anvil/job_outfiles/sock-all_alleleFst_%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

#minInd = 0.3*N or 10, whatever is greater
#minMapQ = 15 bc aligning to a diff species

chrom=NC_068455.1

zcat /home/nhowe/runtiming/anvil/gls/sock-all_${chrom}_minInd0.3.mafs.gz | tail -n +2 -q | cut -f 1,2,3,4 > /home/nhowe/runtiming/anvil/gls/sock-all_${chrom}_minInd0.3.sites
angsd sites index /home/nhowe/runtiming/anvil/gls/sock-all_${chrom}_minInd0.3.sites

angsd -b /home/nhowe/runtiming/anvil/sock-all_EE_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/anvil/gls/sock-all_${chrom}_minInd0.3.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/anvil/diversity/sock-all_${chrom}_EE_minInd0.3 -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 15 -doCounts 1 -minInd 10 -setminDepth 32 -setmaxDepth 640 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/anvil/sock-all_LL_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/anvil/gls/sock-all_${chrom}_minInd0.3.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/anvil/diversity/sock-all_${chrom}_LL_minInd0.3 -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 15 -doCounts 1 -minInd 11 -setminDepth 39 -setmaxDepth 780 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

realSFS /home/nhowe/runtiming/anvil/diversity/sock-all_${chrom}_EE_minInd0.3.saf.idx -fold 1 /home/nhowe/runtiming/anvil/diversity/sock-all_${chrom}_LL_minInd0.3.saf.idx -fold 1 > /home/nhowe/runtiming/anvil/fst/sock-all_${chrom}_EE-LL_minInd0.3.sfs
realSFS fst index /home/nhowe/runtiming/anvil/diversity/sock-all_${chrom}_EE_minInd0.3.saf.idx -fold 1 /home/nhowe/runtiming/anvil/diversity/sock-all_${chrom}_LL_minInd0.3.saf.idx -fold 1 -sfs /home/nhowe/runtiming/anvil/fst/sock-all_${chrom}_EE-LL_minInd0.3.sfs -fstout /home/nhowe/runtiming/anvil/fst/sock-all_${chrom}_EE-LL_minInd0.3.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/anvil/fst/sock-all_${chrom}_EE-LL_minInd0.3.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/anvil/fst/sock-all_${chrom}_EE-LL_minInd0.3.sfs.pbs.fst.txt

