#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --time=0-20:00:00
#SBATCH --job-name=gls_pick-chum
#SBATCH --output=/home/nhowe/runtiming/sock-pick/job_outfiles/pick-chum_gls_minInd0.3_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,standard
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-37%37

module unload bio/angsd/0.933
module load bio/angsd/0.933

JOBS_FILE=/home/nhowe/runtiming/sock-pick/scripts/pick-chum_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	contig=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

angsd -nThreads 8 -b /home/nhowe/runtiming/sock-pick/pick-chum_filtered_bamslist.txt \
	-ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-r ${contig}: \
	-out /home/nhowe/runtiming/sock-pick/gls/pick-chum_${contig}_minInd0.3 \
	-only_proper_pairs 1 -uniqueOnly 1 -remove_bads 1 -trim 0 \
	-C 50 -minMapQ 15 -minQ 20 \
	-doGlf 2 -GL 1 -doMaf 1 -doMajorMinor 1 -doCounts 1 \
	-minMaf 0.05 -SNP_pval 1e-10 -setMinDepth 47 -setMaxDepth 940 -minInd 28
