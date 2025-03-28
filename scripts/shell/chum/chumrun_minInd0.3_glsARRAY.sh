#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=gls0.3chumrun
#SBATCH --output=/home/nhowe/runtiming/chum/job_outfiles/chumrun_minInd0.3_polymorphic_%A-%a.out
#SBATCH --error=/home/nhowe/runtiming/chum/job_outfiles/chumrun_minInd0.3_polymorphic_%A-%a.err
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,standard
#SBATCH --exclude=node[29]
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-37%24

module unload bio/angsd/0.933
module load bio/angsd/0.933

JOBS_FILE=/home/nhowe/runtiming/chum/scripts/chumrun_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	contig=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

# half setminDepth and setmaxDepth
# minInd 0.3*N

angsd -b /home/nhowe/runtiming/chum/chumrun_bamslist.txt -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -r ${contig}: -out /home/nhowe/runtiming/chum/gls/chumrun_${contig}_minInd0.3_minDepthHalf -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -minQ 20 -doCounts 1 -minInd 28 -setminDepth 48 -setmaxDepth 960 -doGlf 2 -GL 1 -doMaf 1 -doMajorMinor 1 -minMaf 0.05 -SNP_pval 1e-10 -doDepth 1 -dumpCounts 3 -only_proper_pairs 1
