#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=pop_coho-chum
#SBATCH --output=/home/nhowe/runtiming/coho/job_outfiles/coho-chum_pop-analyses_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,himem
#SBATCH --exclude=node[29],himem[01-02]
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-37%18

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

JOBS_FILE=/home/nhowe/runtiming/coho/scripts/coho-chum_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	chrom=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

angsd -b /home/nhowe/runtiming/coho/coho-chum_early_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/coho/gls/coho-chum_wholegenome_minInd0.3_minDepthHalf.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/coho/diversity/coho-chum_${chrom}_early_minInd0.3_minDepthHalf -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 15 -doCounts 1 -minInd 13 -setminDepth 22 -setmaxDepth 450 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/coho/coho-chum_late_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/coho/gls/coho-chum_wholegenome_minInd0.3_minDepthHalf.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/coho/diversity/coho-chum_${chrom}_late_minInd0.3_minDepthHalf -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 15 -doCounts 1 -minInd 11 -setminDepth 19 -setmaxDepth 380 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

realSFS /home/nhowe/runtiming/coho/diversity/coho-chum_${chrom}_early_minInd0.3_minDepthHalf.saf.idx -fold 1 /home/nhowe/runtiming/coho/diversity/coho-chum_${chrom}_late_minInd0.3_minDepthHalf.saf.idx -fold 1 > /home/nhowe/runtiming/coho/fst/coho-chum_${chrom}_early-late_minInd0.3_minDepthHalf.sfs
realSFS fst index /home/nhowe/runtiming/coho/diversity/coho-chum_${chrom}_early_minInd0.3_minDepthHalf.saf.idx -fold 1 /home/nhowe/runtiming/coho/diversity/coho-chum_${chrom}_late_minInd0.3_minDepthHalf.saf.idx -fold 1 -sfs /home/nhowe/runtiming/coho/fst/coho-chum_${chrom}_early-late_minInd0.3_minDepthHalf.sfs -fstout /home/nhowe/runtiming/coho/fst/coho-chum_${chrom}_early-late_minInd0.3_minDepthHalf.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/coho/fst/coho-chum_${chrom}_early-late_minInd0.3_minDepthHalf.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/coho/fst/coho-chum_${chrom}_early-late_minInd0.3_minDepthHalf.sfs.pbs.fst.txt

