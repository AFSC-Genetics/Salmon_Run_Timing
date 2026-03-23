#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --job-name=fst_pick-chum
#SBATCH --output=/home/nhowe/runtiming/sock-pick/job_outfiles/pick-chum_fst_analyses_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --partition=standard,medmem
#SBATCH --exclude=node[30]
#SBATCH --array=1-37%24

module purge
module load bio/angsd/0.933

JOBS_FILE=/home/nhowe/runtiming/sock-pick/scripts/pick-chum_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	chrom=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

angsd -nThreads 8 -b /home/nhowe/runtiming/sock-pick/pick-chum_early_bams.txt \
        -r ${chrom}: \
        -sites /home/nhowe/runtiming/sock-pick/gls/pick-chum_wholegenome_minInd0.3.sites \
        -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
        -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-out /home/nhowe/runtiming/sock-pick/diversity/pick-chum_${chrom}_early_minInd0.3 \
        -only_proper_pairs 1 -uniqueOnly 1 -remove_bads 1 -trim 0 \
        -C 50 -minMapQ 15 \
        -doCounts 1 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 \
	-setminDepth 24 -setmaxDepth 480 -minInd 14

angsd -nThreads 8 -b /home/nhowe/runtiming/sock-pick/pick-chum_late_bams.txt \
	-r ${chrom}: \
	-sites /home/nhowe/runtiming/sock-pick/gls/pick-chum_wholegenome_minInd0.3.sites \
	-ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-out /home/nhowe/runtiming/sock-pick/diversity/pick-chum_${chrom}_late_minInd0.3 \
	-only_proper_pairs 1 -uniqueOnly 1 -remove_bads 1 -trim 0 \
	-C 50 -minMapQ 15 \
	-doCounts 1 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 \
	-setminDepth 23 -setmaxDepth 460 -minInd 13

realSFS /home/nhowe/runtiming/sock-pick/diversity/pick-chum_${chrom}_early_minInd0.3.saf.idx -fold 1 \
	/home/nhowe/runtiming/sock-pick/diversity/pick-chum_${chrom}_late_minInd0.3.saf.idx -fold 1 > /home/nhowe/runtiming/sock-pick/fst/pick-chum_${chrom}_early-late_minInd0.3.sfs

realSFS fst index /home/nhowe/runtiming/sock-pick/diversity/pick-chum_${chrom}_early_minInd0.3.saf.idx -fold 1 \
	/home/nhowe/runtiming/sock-pick/diversity/pick-chum_${chrom}_late_minInd0.3.saf.idx -fold 1 \
	-sfs /home/nhowe/runtiming/sock-pick/fst/pick-chum_${chrom}_early-late_minInd0.3.sfs \
	-fstout /home/nhowe/runtiming/sock-pick/fst/pick-chum_${chrom}_early-late_minInd0.3.sfs.pbs -whichFst 1

realSFS fst stats2 /home/nhowe/runtiming/sock-pick/fst/pick-chum_${chrom}_early-late_minInd0.3.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/sock-pick/fst/pick-chum_${chrom}_early-late_minInd0.3.sfs.pbs.fst.txt

