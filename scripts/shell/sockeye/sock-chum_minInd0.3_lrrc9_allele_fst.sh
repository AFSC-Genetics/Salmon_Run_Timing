#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=sock-chum_lrrc9_fst
#SBATCH --output=/home/nhowe/runtiming/sockeye/job_outfiles/sock-chum_alleleFst_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

#minInd = 0.3*N or 10, whatever is greater
#minMapQ = 15 bc aligning to a diff species

chrom=NC_068455.1

angsd -b /home/nhowe/runtiming/sockeye/sock-chum_EE_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/sockeye/gls/sock-chum_wholegenome_minInd0.3.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/sockeye/diversity/sock-chum_${chrom}_EE_minInd0.3 -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -minInd 10 -setminDepth 32 -setmaxDepth 640 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/sockeye/sock-chum_LL_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/sockeye/gls/sock-chum_wholegenome_minInd0.3.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/sockeye/diversity/sock-chum_${chrom}_LL_minInd0.3 -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -minInd 10 -setminDepth 14 -setmaxDepth 280 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

realSFS /home/nhowe/runtiming/sockeye/diversity/sock-chum_${chrom}_EE_minInd0.3.saf.idx -fold 1 /home/nhowe/runtiming/sockeye/diversity/sock-chum_${chrom}_LL_minInd0.3.saf.idx -fold 1 > /home/nhowe/runtiming/sockeye/fst/sock-chum_${chrom}_EE-LL_minInd0.3.sfs
realSFS fst index /home/nhowe/runtiming/sockeye/diversity/sock-chum_${chrom}_EE_minInd0.3.saf.idx -fold 1 /home/nhowe/runtiming/sockeye/diversity/sock-chum_${chrom}_LL_minInd0.3.saf.idx -fold 1 -sfs /home/nhowe/runtiming/sockeye/fst/sock-chum_${chrom}_EE-LL_minInd0.3.sfs -fstout /home/nhowe/runtiming/sockeye/fst/sock-chum_${chrom}_EE-LL_minInd0.3.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/sockeye/fst/sock-chum_${chrom}_EE-LL_minInd0.3.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/sockeye/fst/sock-chum_${chrom}_EE-LL_minInd0.3.sfs.pbs.fst.txt

