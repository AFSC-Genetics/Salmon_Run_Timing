#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=eucl_er1_fst
#SBATCH --output=/home/nhowe/runtiming/anvil/job_outfiles/euclide_pop-analyses_%A.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --partition=medmem,standard

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

chrom="NC_068428.1"

angsd -b /home/nhowe/runtiming/anvil/euclide_EE_er1_bams.txt \
	-r ${chrom}: \
	-sites /home/nhowe/runtiming/anvil/gls/euclide_wholegenome_minInd0.3.sites \
	-ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-out /home/nhowe/runtiming/anvil/diversity/euclide_${chrom}_EE_minInd0.3 \
	-nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 \
	-C 50 -minMapQ 20 -doCounts 1 \
	-GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1 \
	-minInd 8 -setminDepth 13 -setmaxDepth 260

angsd -b /home/nhowe/runtiming/anvil/euclide_LL_er1_bams.txt \
        -r ${chrom}: \
        -sites /home/nhowe/runtiming/anvil/gls/euclide_wholegenome_minInd0.3.sites \
        -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-out /home/nhowe/runtiming/anvil/diversity/euclide_${chrom}_LL_minInd0.3 \
	-nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 \
	-C 50 -minMapQ 20 -doCounts 1 \
	-GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1 \
	-minInd 8 -setminDepth 21 -setmaxDepth 420

realSFS /home/nhowe/runtiming/anvil/diversity/euclide_${chrom}_EE_minInd0.3.saf.idx -fold 1 /home/nhowe/runtiming/anvil/diversity/euclide_${chrom}_LL_minInd0.3.saf.idx -fold 1 > /home/nhowe/runtiming/anvil/fst/euclide_${chrom}_EE-LL_minInd0.3.sfs
realSFS fst index /home/nhowe/runtiming/anvil/diversity/euclide_${chrom}_EE_minInd0.3.saf.idx -fold 1 /home/nhowe/runtiming/anvil/diversity/euclide_${chrom}_LL_minInd0.3.saf.idx -fold 1 -sfs /home/nhowe/runtiming/anvil/fst/euclide_${chrom}_EE-LL_minInd0.3.sfs -fstout /home/nhowe/runtiming/anvil/fst/euclide_${chrom}_EE-LL_minInd0.3.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/anvil/fst/euclide_${chrom}_EE-LL_minInd0.3.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/anvil/fst/euclide_${chrom}_EE-LL_minInd0.3.sfs.pbs.fst.txt

