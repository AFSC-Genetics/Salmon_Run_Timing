#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=counts_euclide
#SBATCH --output=/home/nhowe/runtiming/anvil/job_outfiles/euclide_wholegenome_counts_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=standard,medmem
#SBATCH --exclude=node[01-05,29]
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-37%20

module purge
module load bio/angsd/0.933

# Change paths/filenames as needed:
PROJDIR="/home/nhowe/runtiming/anvil" # ---> change !!!!!!!!!!!!!!!!!!
PREFIX="euclide"
SUFFIX="minInd0.3"
REFGENOME="/home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna"

####################################################################################

JOBS_FILE=${PROJDIR}/scripts/${PREFIX}_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	chrom=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

# minMapQ = 15

angsd -b ${PROJDIR}/${PREFIX}_early_bams.txt -r ${chrom}: \
	-sites ${PROJDIR}/gls/${PREFIX}_wholegenome_${SUFFIX}.sites \
	-ref ${REFGENOME} -anc ${REFGENOME} \
	-out ${PROJDIR}/diversity/${PREFIX}_${chrom}_early_${SUFFIX}_MM3 \
	-nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -only_proper_pairs 1 \
	-minMapQ 15 -minInd 8 -setminDepth 27 -setmaxDepth 540 \
	-GL 1 -doMaf 1 -doMajorMinor 3 -doCounts 1 -dumpCounts 3

angsd -b ${PROJDIR}/${PREFIX}_late_bams.txt -r ${chrom}: \
	-sites ${PROJDIR}/gls/${PREFIX}_wholegenome_${SUFFIX}.sites \
	-ref ${REFGENOME} -anc ${REFGENOME} \
	-out ${PROJDIR}/diversity/${PREFIX}_${chrom}_late_${SUFFIX}_MM3 \
	-nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -only_proper_pairs 1 \
	-minMapQ 15 -minInd 8 -setminDepth 27 -setmaxDepth 540 \
	-GL 1 -doMaf 1 -doMajorMinor 3 -doCounts 1 -dumpCounts 3

