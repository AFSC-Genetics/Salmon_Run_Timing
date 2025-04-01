#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=pop_chumrun
#SBATCH --output=/home/nhowe/runtiming/chum/job_outfiles/chumrun_esrb_fst_%A.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

# minInd 0.3*N
# minDepth = 0.5*N
# minMap Q = 20

chrom=NC_068449.1

angsd -b /home/nhowe/runtiming/chum/chumrun_EE_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/chum/gls/chumrun_wholegenome_minInd0.3_minDepthHalf.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_EE_minInd0.3_minDepthHalf -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -minInd 10 -setminDepth 19 -setmaxDepth 390 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/chum/chumrun_LL_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/chum/gls/chumrun_wholegenome_minInd0.3_minDepthHalf.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_LL_minInd0.3_minDepthHalf -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -minInd 16 -setminDepth 23 -setmaxDepth 470 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

realSFS /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_EE_minInd0.3_minDepthHalf.saf.idx -fold 1 /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_LL_minInd0.3_minDepthHalf.saf.idx -fold 1 > /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_EE-LL_minInd0.3_minDepthHalf.sfs
realSFS fst index /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_EE_minInd0.3_minDepthHalf.saf.idx -fold 1 /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_LL_minInd0.3_minDepthHalf.saf.idx -fold 1 -sfs /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_EE-LL_minInd0.3_minDepthHalf.sfs -fstout /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_EE-LL_minInd0.3_minDepthHalf.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt

