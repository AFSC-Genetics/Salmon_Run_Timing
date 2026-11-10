#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=fis_AukePink
#SBATCH --output=/share2/nhowe/runtiming/pink/home/pbarry/AukePink_lcWGS/job_outfiles/AukePink_minInd0.3_minMaf_fis_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=standard,medmem
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=10

module purge
module load bio/angsd/0.933

REFGENOME='/home/nhowe/reference_genomes/pink_even/GCF_021184085.1_OgorEven_v1.0_genomic.fna'

JOBS_FILE=/share2/nhowe/runtiming/pink/home/pbarry/AukePink_lcWGS/scripts/AukePink_chromARRAY.txt
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

angsd -nThreads 10 -b /share2/nhowe/runtiming/pink/home/pbarry/AukePink_lcWGS/AukePink_EarlyEven_bams_NH.txt \
	-r ${chrom}: \
	-sites /share2/nhowe/runtiming/pink/home/pbarry/AukePink_lcWGS/gls/AukePink_wholegenome_polymorphic.sites \
	-ref ${REFGENOME} \
	-anc ${REFGENOME} \
	-out /share2/nhowe/runtiming/fis/AukePink_${chrom}_early_minInd0.3_HWE_minMaf \
	-uniqueOnly 1 -remove_bads 1 -trim 0 -only_proper_pairs 1 \
	-C 50 -minMapQ 20 \
	-GL 1 -dumpCounts 3 -doMaf 1 -doMajorMinor 3 -doHWE 1 -doCounts 1 \
	-minInd 22 -setminDepth 73 -setmaxDepth 1460 \
	-minMaf 0.05 -SNP_pval 1e-10

angsd -nThreads 10 -b /share2/nhowe/runtiming/pink/home/pbarry/AukePink_lcWGS/AukePink_LateEven_bams_NH.txt \
        -r ${chrom}: \
        -sites /share2/nhowe/runtiming/pink/home/pbarry/AukePink_lcWGS/gls/AukePink_wholegenome_polymorphic.sites \
        -ref ${REFGENOME} \
        -anc ${REFGENOME} \
	-out /share2/nhowe/runtiming/fis/AukePink_${chrom}_late_minInd0.3_HWE_minMaf \
        -uniqueOnly 1 -remove_bads 1 -trim 0 -only_proper_pairs 1 \
        -C 50 -minMapQ 20 \
        -GL 1 -dumpCounts 3 -doMaf 1 -doMajorMinor 3 -doHWE 1 -doCounts 1 \
	-minInd 22 -setminDepth 74 -setmaxDepth 1480 \
        -minMaf 0.05 -SNP_pval 1e-10


