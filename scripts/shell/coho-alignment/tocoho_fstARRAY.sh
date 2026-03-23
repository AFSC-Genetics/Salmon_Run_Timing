#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=pop_tocoho
#SBATCH --output=/home/nhowe/runtiming/tocoho/job_outfiles/tocoho_pop-analyses_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=medmem,standard,himem
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-30%15

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

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

angsd -b /home/nhowe/runtiming/tocoho/tocoho_early_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_polymorphic.sites -ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -anc /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -out /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_early_polymorphic -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -setminDepth 45 -setmaxDepth 900.0 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/tocoho/tocoho_late_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/tocoho/gls/tocoho_wholegenome_polymorphic.sites -ref /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -anc /home/nhowe/reference_genomes/cohoV2/GCF_002021735.2_Okis_V2_genomic.fna -out /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_late_polymorphic -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 20 -doCounts 1 -setminDepth 38 -setmaxDepth 760.0 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

realSFS /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_early_polymorphic.saf.idx -fold 1 /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_late_polymorphic.saf.idx -fold 1 > /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_polymorphic.sfs
realSFS fst index /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_early_polymorphic.saf.idx -fold 1 /home/nhowe/runtiming/tocoho/diversity/tocoho_${chrom}_late_polymorphic.saf.idx -fold 1 -sfs /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_polymorphic.sfs -fstout /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_polymorphic.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_polymorphic.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/tocoho/fst/tocoho_${chrom}_early-late_polymorphic.sfs.pbs.fst.txt

