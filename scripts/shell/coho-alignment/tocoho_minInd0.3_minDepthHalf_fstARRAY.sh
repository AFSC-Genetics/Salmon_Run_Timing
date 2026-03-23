#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=fst_tocoho
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/tocoho_minInd0.3_minDepthHalf_pop-analyses_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=standard,medmem
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-30%15

module purge
module load bio/angsd/0.933


mkdir -p /home/nhowe/runtiming/tocoho/diversity
mkdir -p /home/nhowe/runtiming/tocoho/fst

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

angsd -nThreads 10 -b /home/nhowe/runtiming/tocoho/tocoho_early_bams_scratch2.txt \
	-r ${chrom}: \
	-sites /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minInd0.3_minDepthHalf.sites \
	-ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna \
	-anc /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna \
	-out /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_early_minInd0.3_minDepthHalf \
	-uniqueOnly 1 -remove_bads 1 -trim 0 -only_proper_pairs 1 \
	-C 50 -minMapQ 20 \
	-GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -doCounts 1 \
	-minInd 13 -setminDepth 22 -setmaxDepth 450

angsd -nThreads 10 -b /home/nhowe/runtiming/tocoho/tocoho_late_bams_scratch2.txt \
        -r ${chrom}: \
        -sites /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_minInd0.3_minDepthHalf.sites \
        -ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna \
        -anc /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna \
	-out /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_late_minInd0.3_minDepthHalf \
        -uniqueOnly 1 -remove_bads 1 -trim 0 -only_proper_pairs 1 \
        -C 50 -minMapQ 20 \
        -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -doCounts 1 \
	-minInd 10 -setminDepth 16 -setmaxDepth 380

realSFS /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_early_minInd0.3_minDepthHalf.saf.idx /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_late_minInd0.3_minDepthHalf.saf.idx -fold 1 > /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_minInd0.3_minDepthHalf.sfs

realSFS fst index /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_early_minInd0.3_minDepthHalf.saf.idx /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_late_minInd0.3_minDepthHalf.saf.idx -fold 1 \
	-sfs /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_minInd0.3_minDepthHalf.sfs \
	-fstout /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_minInd0.3_minDepthHalf.sfs.pbs -whichFst 1

realSFS fst stats2 /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_minInd0.3_minDepthHalf.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_minInd0.3_minDepthHalf.sfs.pbs.fst.txt

