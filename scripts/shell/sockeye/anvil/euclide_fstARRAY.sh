#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --job-name=pop_euclide
#SBATCH --output=/home/nhowe/runtiming/anvil/job_outfiles/euclide_pop-analyses_%A-%a.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=natasha.howe@noaa.gov
#SBATCH --array=1-37%24

module unload bio/angsd/0.933 bio/ngstools/202202
module load bio/angsd/0.933 bio/ngstools/202202

JOBS_FILE=/home/nhowe/runtiming/anvil/scripts/euclide_angsdARRAY_input.txt
IDS=$(cat ${JOBS_FILE})

for sample_line in ${IDS}
do
	job_index=$(echo ${sample_line} | awk -F ":" '{print $1}')
	chrom=$(echo ${sample_line} | awk -F ":" '{print $2}')
	if [[ ${SLURM_ARRAY_TASK_ID} == ${job_index} ]]; then
		break
	fi
done

angsd -b /home/nhowe/runtiming/anvil/euclide_beach_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/anvil/gls/euclide_wholegenome_minInd0.3.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/anvil/diversity/euclide_${chrom}_beach_minInd0.3 -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 15 -doCounts 1 -minInd 8 -setminDepth 27 -setmaxDepth 540.0 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

angsd -b /home/nhowe/runtiming/anvil/euclide_creek_bams.txt -r ${chrom}: -sites /home/nhowe/runtiming/anvil/gls/euclide_wholegenome_minInd0.3.sites -ref /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -anc /home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna -out /home/nhowe/runtiming/anvil/diversity/euclide_${chrom}_creek_minInd0.3 -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 0 -C 50 -minMapQ 15 -doCounts 1 -minInd 8 -setminDepth 27 -setmaxDepth 540.0 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doSaf 1 -only_proper_pairs 1

realSFS /home/nhowe/runtiming/anvil/diversity/euclide_${chrom}_beach_minInd0.3.saf.idx -fold 1 /home/nhowe/runtiming/anvil/diversity/euclide_${chrom}_creek_minInd0.3.saf.idx -fold 1 > /home/nhowe/runtiming/anvil/fst/euclide_${chrom}_beach-creek_minInd0.3.sfs
realSFS fst index /home/nhowe/runtiming/anvil/diversity/euclide_${chrom}_beach_minInd0.3.saf.idx -fold 1 /home/nhowe/runtiming/anvil/diversity/euclide_${chrom}_creek_minInd0.3.saf.idx -fold 1 -sfs /home/nhowe/runtiming/anvil/fst/euclide_${chrom}_beach-creek_minInd0.3.sfs -fstout /home/nhowe/runtiming/anvil/fst/euclide_${chrom}_beach-creek_minInd0.3.sfs.pbs -whichFst 1
realSFS fst stats2 /home/nhowe/runtiming/anvil/fst/euclide_${chrom}_beach-creek_minInd0.3.sfs.pbs.fst.idx -win 1 -step 1 > /home/nhowe/runtiming/anvil/fst/euclide_${chrom}_beach-creek_minInd0.3.sfs.pbs.fst.txt

