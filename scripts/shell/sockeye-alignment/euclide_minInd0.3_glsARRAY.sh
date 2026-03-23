#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=0-20:00:00
#SBATCH --job-name=gls_sockeye
#SBATCH --output=/home/nhowe/runtiming/sock-sock/job_outfiles/euclide_minInd0.3_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-29%13

module purge
module load bio/angsd/0.933

# SOCKEYE V1
REFGENOME="/home/nhowe/reference_genomes/sockeyeV1/GCF_006149115.2_Oner_1.1_genomic.fna"
OUTDIR="/scratch2/nhowe/runtiming/sock-sock/gls"
echo "outdir: "${OUTDIR}

JOBS_FILE=/home/nhowe/maturity/sockeye/scripts/sockeye_mat_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	contig=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

angsd -b /home/nhowe/runtiming/sock-sock/euclide_filtered_bamslist.txt -ref ${REFGENOME} \
	-r ${contig}: \
	-out ${OUTDIR}/euclide_${contig}_minInd0.3 \
	-nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -only_proper_pairs 1 \
	-C 50 -minMapQ 20 -minQ 20 \
	-doCounts 1 -dumpCounts 3 \
	-minInd 16 -setminDepth 54 -setmaxDepth 1080 \
	-doGlf 2 -GL 1 -doMaf 1 -doMajorMinor 1 \
	-minMaf 0.05 -SNP_pval 1e-10 -doDepth 1
