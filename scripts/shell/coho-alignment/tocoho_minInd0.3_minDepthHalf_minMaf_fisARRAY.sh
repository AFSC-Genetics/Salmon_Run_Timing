#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=fis_tocoho
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/tocoho_minInd0.3_minDepthHalf_minMaf_fis_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=standard,medmem
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=7,21

### --- CHANGED OUTPUT LOCATION TO SHARE2 BC OF SPACE ISSUES --- ####

### also, only need this for the chromosomes of interest now #####

module purge
module load bio/angsd/0.933

JOBS_FILE=/home/nhowe/runtiming/tocoho/scripts/tocoho_chromARRAY.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	chrom=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

#minInd 30% minDepth Half
# added minmaf and snp_pval back because this is not a SFS, and we need to remove monomorphic sites

angsd -nThreads 10 -b /home/nhowe/runtiming/tocoho/tocoho_early_bams.txt \
	-r ${chrom}: \
	-sites /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minInd0.3_minDepthHalf.sites \
	-ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna \
	-anc /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna \
	-out /share2/nhowe/runtiming/fis/tocoho_${chrom}_early_minInd0.3_minDepthHalf_HWE_minMaf \
	-uniqueOnly 1 -remove_bads 1 -trim 0 -only_proper_pairs 1 \
	-C 50 -minMapQ 20 \
	-GL 1 -dumpCounts 3 -doMaf 1 -doMajorMinor 3 -doHWE 1 -doCounts 1 \
	-minInd 13 -setminDepth 22 -setmaxDepth 450 \
	-minMaf 0.05 -SNP_pval 1e-10

angsd -nThreads 10 -b /home/nhowe/runtiming/tocoho/tocoho_late_bams.txt \
        -r ${chrom}: \
        -sites /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minInd0.3_minDepthHalf.sites \
        -ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna \
        -anc /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna \
	-out /share2/nhowe/runtiming/fis/tocoho_${chrom}_late_minInd0.3_minDepthHalf_HWE_minMaf \
        -uniqueOnly 1 -remove_bads 1 -trim 0 -only_proper_pairs 1 \
        -C 50 -minMapQ 20 \
        -GL 1 -dumpCounts 3 -doMaf 1 -doMajorMinor 3 -doHWE 1 -doCounts 1 \
	-minInd 10 -setminDepth 16 -setmaxDepth 380 \
        -minMaf 0.05 -SNP_pval 1e-10


