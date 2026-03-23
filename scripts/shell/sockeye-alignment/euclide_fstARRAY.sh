#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=fst_euclide
#SBATCH --partition=medmem,himem,standard
#SBATCH --output=/home/nhowe/runtiming/sock-sock/job_outfiles/euclide_fst_analyses_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-37%12

module purge
module load bio/angsd/0.933

REFGENOME="/home/nhowe/reference_genomes/sockeyeV1/GCF_006149115.2_Oner_1.1_genomic.fna"

JOBS_FILE=/home/nhowe/maturity/sockeye/scripts/sockeye_mat_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	chrom=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

angsd -b /home/nhowe/runtiming/sock-sock/euclide_late_bams.txt \
	-r ${chrom}: \
	-sites /scratch2/nhowe/runtiming/sock-sock/gls/euclide_wholegenome_minInd0.3.sites \
	-ref ${REFGENOME} -anc ${REFGENOME} \
	-out /scratch2/nhowe/runtiming/sock-sock/diversity/euclide_${chrom}_late_minInd0.3 \
	-nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 \
	-doCounts 1 -dumpCounts 3 \
	-minInd 8 -setminDepth 27 -setmaxDepth 540 \
	-GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/sock-sock/euclide_early_bams.txt \
	-r ${chrom}: -sites /scratch2/nhowe/runtiming/sock-sock/gls/euclide_wholegenome_minInd0.3.sites \
	-ref ${REFGENOME} -anc ${REFGENOME} \
	-out /scratch2/nhowe/runtiming/sock-sock/diversity/euclide_${chrom}_early_minInd0.3 \
	-nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 \
	-doCounts 1 -dumpCounts 3 \
	-minInd 8 -setminDepth 27 -setmaxDepth 540 \
	-GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

realSFS /scratch2/nhowe/runtiming/sock-sock/diversity/euclide_${chrom}_late_minInd0.3.saf.idx /scratch2/nhowe/runtiming/sock-sock/diversity/euclide_${chrom}_early_minInd0.3.saf.idx -fold 1 > /scratch2/nhowe/runtiming/sock-sock/fst/euclide_${chrom}_late-early_minInd0.3.sfs

realSFS fst index /scratch2/nhowe/runtiming/sock-sock/diversity/euclide_${chrom}_late_minInd0.3.saf.idx /scratch2/nhowe/runtiming/sock-sock/diversity/euclide_${chrom}_early_minInd0.3.saf.idx -fold 1 \
	-sfs /scratch2/nhowe/runtiming/sock-sock/fst/euclide_${chrom}_late-early_minInd0.3.sfs -whichFst 1 \
	-fstout /scratch2/nhowe/runtiming/sock-sock/fst/euclide_${chrom}_late-early_minInd0.3.sfs.pbs

realSFS fst stats2 /scratch2/nhowe/runtiming/sock-sock/fst/euclide_${chrom}_late-early_minInd0.3.sfs.pbs.fst.idx \
	-win 1 -step 1 > /scratch2/nhowe/runtiming/sock-sock/fst/euclide_${chrom}_late-early_minInd0.3.sfs.pbs.fst.txt

