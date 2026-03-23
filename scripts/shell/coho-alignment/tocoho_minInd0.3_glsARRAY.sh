#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=plm_tocoho
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/tocoho_minInd0.3_minDepthHalf_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,standard
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-30%15

module unload bio/angsd/0.933
module load bio/angsd/0.933


JOBS_FILE=/home/nhowe/runtiming/tocoho/scripts/tocoho_chromARRAY.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	contig=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

# 0.4x cutoff
# minInd = 0.3*N
# minDepth = 0.5*N
# maxDepth = 10*N

angsd -nThreads 10 -b /home/nhowe/runtiming/tocoho/tocoho_filtered_bamslist_scratch.txt \
	-ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna \
	-r ${contig}: \
	-out /home/nhowe/runtiming/tocoho/gls/tocoho_${contig}_minInd0.3_minDepthHalf \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 \
	-C 50 -minMapQ 20 -minQ 20 \
	-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 \
	-minMaf 0.05 -SNP_pval 1e-10 \
	-doCounts 1 -dumpCounts 3 -doDepth 1 \
	-setminDepth 41 -setmaxDepth 830 -minInd 24
