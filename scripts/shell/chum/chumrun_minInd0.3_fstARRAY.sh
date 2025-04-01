#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=Ind0.3Fst_chumrun
#SBATCH --output=/home/nhowe/runtiming/chum/job_outfiles/chumrun_minInd0.3_fst_fall-summer_%A-%a.out
#SBATCH --error=/home/nhowe/runtiming/chum/job_outfiles/chumrun_minInd0.3_fst_fall-summer_%A-%a.err
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-37%24

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

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

angsd -b /home/nhowe/runtiming/chum/chumrun_fall_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/chum/gls/chumrun_wholegenome_minInd0.3_minDepthHalf.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_fall_minInd0.3_minDepthHalf -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -minInd 14 -setminDepth 24 -setmaxDepth 480 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/chum/chumrun_summer_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/chum/gls/chumrun_wholegenome_minInd0.3_minDepthHalf.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_summer_minInd0.3_minDepthHalf -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -minInd 14 -setminDepth 24 -setmaxDepth 480 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

realSFS /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_fall_minInd0.3_minDepthHalf.saf.idx -fold 1 /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_summer_minInd0.3_minDepthHalf.saf.idx -fold 1 > /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_fall-summer_minInd0.3_minDepthHalf.sfs
realSFS fst index /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_fall_minInd0.3_minDepthHalf.saf.idx -fold 1 /home/nhowe/runtiming/chum/diversity/chumrun_${chrom}_summer_minInd0.3_minDepthHalf.saf.idx -fold 1 -sfs /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_fall-summer_minInd0.3_minDepthHalf.sfs -fstout /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_fall-summer_minInd0.3_minDepthHalf.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_fall-summer_minInd0.3_minDepthHalf.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/chum/fst/chumrun_${chrom}_fall-summer_minInd0.3_minDepthHalf.sfs.pbs.fst.txt

