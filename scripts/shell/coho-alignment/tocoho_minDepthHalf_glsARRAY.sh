#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=plm_tocoho
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/tocoho_minDepthHalf_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,himem,standard
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-30%20

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

# min depth half, max depth 10* sample size

angsd -b /home/nhowe/runtiming/tocoho/tocoho_filtered_bamslist1.txt -ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -r ${contig}: -out /home/nhowe/runtiming/tocoho/gls/tocoho_${contig}_minDepthHalf -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -minQ 20 -doCounts 1 -setminDepth 41 -setmaxDepth 830 -doGlf 2 -GL 1 -doMaf 1 -doMajorMinor 1 -minMaf 0.05 -SNP_pval 1e-10 -doDepth 1 -dumpCounts 3 -only_proper_pairs 1
