#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=fis_chumrun
#SBATCH --output=/home/nhowe/runtiming/chum/job_outfiles/chumrun_minInd0.3_minDepthHalf_fis_minMaf_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=standard,medmem
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=8,29,35

module purge
module load bio/angsd/0.933

JOBS_FILE=/home/nhowe/runtiming/chum/scripts/chumrun_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	chrom=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

# setminDepth = ~50% sample size
# setmaxDepth = 10 x sample size
# minInd = 0.3*samplesize
#minMapQ = 20
# added  -minMaf 0.05 -SNP_pval 1e-10 bc not using SFS and instead need only biallelic sites for FIS calculated independently per run timing


angsd -b /home/nhowe/runtiming/chum/chumrun_fall_bams.txt -r ${chrom}: \
	-sites /home/nhowe/runtiming/chum/gls/chumrun_wholegenome_minInd0.3_minDepthHalf.sites \
	-ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-out /share2/nhowe/runtiming/fis/chumrun_${chrom}_fall_minInd0.3_minDepthHalf_HWE_minMaf \
	-nThreads 10  -only_proper_pairs 1 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 \
	-doCounts 1 -minInd 14 -setminDepth 24 -setmaxDepth 480 -GL 1 -doMaf 1 -doMajorMinor 3 \
	-doHWE 1 -minMaf 0.05 -SNP_pval 1e-10

angsd -b /home/nhowe/runtiming/chum/chumrun_summer_bams.txt -r ${chrom}: \
	-sites /home/nhowe/runtiming/chum/gls/chumrun_wholegenome_minInd0.3_minDepthHalf.sites \
	-ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna \
	-out /share2/nhowe/runtiming/fis/chumrun_${chrom}_summer_minInd0.3_minDepthHalf_HWE_minMaf \
	-nThreads 10 -only_proper_pairs 1 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 \
	-doCounts 1 -minInd 14 -setminDepth 24 -setmaxDepth 480 -GL 1 -doMaf 1 -doMajorMinor 3 \
	-doHWE 1 -minMaf 0.05 -SNP_pval 1e-10

