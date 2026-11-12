#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:00
#SBATCH --job-name=plm_AukePink
#SBATCH --output=/home/pbarry/AukePink_lcWGS/job_outfiles/AukePink_polymorphic_%A-%a.out
#SBATCH --array=1-27%24

module unload bio/angsd/0.933
module load bio/angsd/0.933

JOBS_FILE=/home/pbarry/AukePink_lcWGS/scripts/AukePink_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	contig=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

angsd -b /home/pbarry/AukePink_lcWGS/AukePink_filtered_bamslist.txt -ref /home/pbarry/ref_genomes/GCF_021184085.1_OgorEven_v1.0_genomic.fna -r ${contig}: -out /home/pbarry/AukePink_lcWGS/gls/AukePink_${contig}_polymorphic -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 15 -doCounts 1 -setminDepth 276 -setmaxDepth 5520.0 -doGlf 2 -GL 1 -doMaf 1 -doMajorMinor 1 -minMaf 0.05 -SNP_pval 1e-10 -doDepth 1 -dumpCounts 3 -only_proper_pairs 1