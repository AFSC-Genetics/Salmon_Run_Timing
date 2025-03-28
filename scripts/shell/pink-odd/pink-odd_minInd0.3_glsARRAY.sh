#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=gls_pink-odd
#SBATCH --output=/home/nhowe/runtiming/pinkOdd/job_outfiles/pink-odd_gls_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-37%37

module unload bio/angsd/0.933
module load bio/angsd/0.933

JOBS_FILE=/home/nhowe/runtiming/pinkOdd/scripts/pink-odd_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	contig=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

# minInd 0.3

angsd -b /home/nhowe/runtiming/pinkOdd/pink-odd_filtered_bamslist.txt -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -r ${contig}: -out /home/nhowe/runtiming/pinkOdd/gls/pink-odd_${contig}_minInd0.3 -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 15 -minQ 20 -doCounts 1 -minInd 40 -setminDepth 135 -setmaxDepth 2700 -doGlf 2 -GL 1 -doMaf 1 -doMajorMinor 1 -minMaf 0.05 -SNP_pval 1e-10 -doDepth 1 -dumpCounts 3 -only_proper_pairs 1
