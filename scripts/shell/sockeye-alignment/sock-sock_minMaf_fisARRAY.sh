#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=fis_sock
#SBATCH --partition=medmem,standard
#SBATCH --output=/home/nhowe/runtiming/sock-sock/job_outfiles/sock-sock_minInd0.3_minMaf_fis_analyses_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=12,13

### --- MOVED OUTPUTS TO SHARE2 BC OF SPACE LIMITS --- ###

module purge
module load bio/angsd/0.933

REFGENOME="/home/nhowe/reference_genomes/sockeyeV1/GCF_006149115.2_Oner_1.1_genomic.fna"

mkdir -p /home/nhowe/runtiming/sock-sock/diversity

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

# added  -minMaf 0.05 -SNP_pval 1e-10 bc not running SFS and therefore need only biallelic sites (drop monomorphic w/in a run)

angsd -b /home/nhowe/runtiming/sock-sock/euclide_late_bams_anvilpath.txt \
	-r ${chrom}: \
	-sites /home/nhowe/runtiming/sock-sock/gls/euclide_wholegenome_minInd0.3.sites \
	-ref ${REFGENOME} -anc ${REFGENOME} \
	-out /share2/nhowe/runtiming/fis/sock-sock_${chrom}_late_minInd0.3_HWE_minMaf \
	-nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 \
	-doCounts 1 -dumpCounts 3 -doHWE 1 \
	-minInd 8 -setminDepth 27 -setmaxDepth 540 \
	-GL 1 -doMaf 1 -doMajorMinor 3 -only_proper_pairs 1 \
	-minMaf 0.05 -SNP_pval 1e-10

angsd -b /home/nhowe/runtiming/sock-sock/euclide_early_bams_anvilpath.txt \
	-r ${chrom}: -sites /home/nhowe/runtiming/sock-sock/gls/euclide_wholegenome_minInd0.3.sites \
	-ref ${REFGENOME} -anc ${REFGENOME} \
	-out /share2/nhowe/runtiming/fis/sock-sock_${chrom}_early_minInd0.3_HWE_minMaf \
	-nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 \
	-doCounts 1 -dumpCounts 3 -doHWE 1 \
	-minInd 8 -setminDepth 27 -setmaxDepth 540 \
	-GL 1 -doMaf 1 -doMajorMinor 3 -only_proper_pairs 1 \
	 -minMaf 0.05 -SNP_pval 1e-10
