---
title: "Runtiming - Four Salmon Species"
author: "Patrick Barry & Natasha Howe"
date: "April 01, 2025"
output: 
  html_document:
    keep_md: true
    wrap: 72
---



## Introduction

Investigating run timing differentiation between early and late returning salmon from multiple regions of Alaska. All species were aligned to chum to compare the specified regions on the genome that are differentiated across species.

## Sample Collection

Data from Auke Creek pink salmon (even and odd lineage), Wood River sockeye salmon, Yukon River chum salmon, and Skilak River coho salmon. Sample collection data was compiled from the raw data in each four species manually with three columns: sampleID, Runtime, and Species


``` r
# All samples and sequencing information can be found at:
read.csv("./data/raw/fourspecies_runtiming_metadata.csv")
```

### Site Map

@PatrickBarry-NOAA


``` r
FigureS1_SiteMap.R
```

## Low-coverage Whole Genome Sequencing Pipeline

Sequencing facilities differed across samples and had vastly different coverages. Refer to methods for species-specific sequencing effort details.

### Step 0: Setup

PREFIX below can be substituted depending on the species comparison being looked at: For example:

Prefix Options:

1.  **Chum (Yukon):** "chumrun"
2.  **Coho (Skilak Lake):** "coho-chum"
3.  **Sockeye**:
    i.  **Teal (Early Creek) & Whitefish (Late Creek)**: "sock-chum"
    ii. **Teal (Early Creek) & Anvil (Late Beach)**: "euclide"
    iii. **Teal (Early Creek) & Whitefish (Late Creek) & Anvil (Late Beach)**: "sock-all"
4.  **Pink:**
    i.  **Even**: "pink-chum"
    ii. **Odd**: "pink-odd

We first downloaded the reference genome onto our Slurm Manager so that we can align our data. We are using the *chum salmon* reference genome v2.

We listed chromosomes by using a few commands in the terminal.


``` bash
grep '^>NC' ./home/nhowe/reference_genomes/chum/GCF_023373465.1_Oket_V2_genomic.fna |  awk '{print $1}' | sed 's/>//' > chromosomes_all.txt
```

Chum have 37 chromosomes + 1 mtDNA.

Remove the mitochondrial genome from our list of chromosomes to analyze.


``` r
sed '/NC_017838.1/d' chromosomes_all.txt > chromosomes.txt
```

Get sizes of chromosomes.


``` bash
 cut -f1,2 GCF_023373465.1_Oket_V2_genomic.fna.fai | grep '^NC' > chrom_sizes.txt
```

Two shell scripts index the reference genome and can be run with the sbatch command which will have the SLURM manager schedule them.


``` bash
sbatch GCF_023373465.1_Oket_V2_genomic_bwa-index_script.sh
sbatch GCF_023373465.1_Oket_V2_genomic_fai_script.sh
```

### Step 1: Quality Control

This step is a quality control of the sequencing data (raw fastqs). It runs FastQC [@Andrews2010] on the raw sequence data. FastQC is a java based application that runs a set of quality control measures.

-   Import data from BAM, SAM or FastQ files (any variant)
-   Provides a quick overview indicating in which areas there may be problems
-   Summary graphs and tables to quickly assess data quality
-   Export of results to an HTML based permanent report
-   Offline operation to allow automated generation of reports without running the interactive application

Run two shell scripts [PREFIX]-raw_fastqcARRAY.sh and [PREFIX]-raw_multiqcSLURM.sh in succession.


``` bash
sbatch [PREFIX]-raw_fastqcARRAY.sh
sbatch [PREFIX]-raw_multiqcSLURM.sh
```

### Step 2: Trim

Trim adapters from raw fastqs with the program TRIMMOMATIC [@Bolger2014] and quality-check the trimmed fastqs with multiQC [@Ewels2016]. Run three shell scripts [PREFIX]\_trimARRAY.sh, [PREFIX]-trim_fastqcARRAY.sh, and [PREFIX]-trim_multiqcSLURM.sh which should all run in succession.


``` bash
# trim adapters
sbatch [PREFIX]_trimARRAY.sh

# multiqc report
sbatch [PREFIX]-trim_fastqcARRAY.sh
sbatch [PREFIX]-trim_multiqcSLURM.sh
```

Download the multiQC.html file and evaluate the data quality before proceeding.

### Step 3: Align

Align reads to the reference genome that we downloaded with BWA and runs the aligned reads through SAMTOOLS [@Li2009]: 'fixmate' cleans up the read pairings and flags from BWA; a pair of 'view' statements converts the .sam file to a .bam file and filters the .bam file for non-unique and poor quality mappings; and 'sort' sorts the read pairings by coordinate (instead of read name). After a .bam file is built, duplicate reads are removed and (if the data is PE) overlapping reads are clipped to generate the final .bam.


``` bash
# align to genome
sbatch [PREFIX]_alignARRAY.sh

# calculate average depths
sbatch [PREFIX]_depthsARRAY.sh
```

There is an output file [PREFIX]\_depths.csv which we can import and run through the downsampling script to see if there is uneven distribution of depths across early and late individuals. We generally used 0.5x as a depth cutoff. Two chum samples were close to the 0.5x cutoff, above or rounding up to 0.4x; therefore, we didn't use a cutoff for chum.


``` r
./species-specific/[PREFIX]_depth_and_downsampling.R
```

#### Mapped Reads

Calculate the total number of mapped reads from the bamslists with samtools.


``` bash
sbatch [PREFIX]_mapped_reads.sh
```

Now calculate the average and standard deviation for the mapped reads after downloading.


``` r
mean(read.csv("./results/depth/[PREFIX]_mapped_reads.csv", header = F)$V1)
sd(read.csv("./results/depth/[PREFIX]_mapped_reads.csv", header = F)$V1)
```

### Step 4: Genotype Likelihoods

Calculates genotype likelihoods for putatively polymorphic sites (SNPs). This is done on a chrom-by-chrom level using an array. These script uses [PREFIX]\_angsdARRAY_input.txt (array file with chromosomes) and a list of the bam files to analyze ([PREFIX]\_filtered_bamslist.txt).


``` bash
sbatch [PREFIX]_minInd0.3_[minDepthHalf_]glsARRAY.sh 
```

Filters slightly differed depending on the species alignment and average coverage.

-   *minMapQ* [15 or 20]: the minimum map quality is 15, which was lowered from 20 for all except chum to allow for additional misalignments since non-chum samples are aligned to a different species.
-   *minQ* [20]: the minimum quality of nucleotide call is 20.
-   *C* [50]: use an adjustment of the Map Quality for excessive mismatches of 50. This is the suggested value for BWA and leads to fewer false positive variant calls.
-   *minDepth* [N or 0.5\*N]: the minimum depth is set the number of individuals if average depth \> 1 and to half the number of individuals if avg depth was \< 1. Discard site if total sequencing depth (all individuals added together) is below minDepth.
-   *setMaxDepth* [20\*minDepth]: maximum depth is set to the number of individuals times 10 or 20. Discard site if total sequencing depth (all individuals added together) is above the number of individuals multiplied by 10 or 20 depending on the average depth of coverage.
-   *SNP_pval* [1e-10]: only work with sites with a p-value from the likelihood test of less than 1e-10.
-   *Gl* [1]: use SAMtools genotype likelihood calling algorithm. One assumption of the SAMtools algorithm is that it assumes that errors are not independent, but that once a first error occurs at a certain site in an individual, a second error is more likely to occur at the same site.
-   *minMaf* [0.05]: only work with sites with a maf above 0.05.
-   *minInd* [0.3\*N]: Set to 30% of sample size to include a proportion of individuals.
-   Additional Filters: *only_proper_pairs* [1]; *remove_bads* [1]; *uniqueOnly* [1]; *trim* [0]; *doMajorMinor* [1]

### Step 5: Collate

The genotype likelihood files from above are split by chromosome. Create whole genome beagle and maf files by concatenating the chromosome beagle and maf files. Also create a *sites* file that will be used as input for Fst calculations under the -sites flag.


``` bash
sbatch [PREFIX]_minInd0.3_concatenate_mafs.sh
sbatch [PREFIX]_minInd0.3_concatenate_beagles.sh
```

### Step 6: Population Analyses - FST

Created bamslists for the early and late run timing phenotypes (subset the filtered bamslist of all samples for that species by the run timing phenotype).

First, angsd calls GLs for each subgroup, but it uses the *sites* flag from the sites file created in Step 5, so the SNPs being called will only be from the polymorphic files from Step 4. Therefore, some of the filtering steps have been removed, such as: *SNP_pval* and *minMaf*. Filters were adjusted for the sample size of the subset group: *setminDepth* and set*maxDepth*, and *minInd.*


``` bash
sbatch [PREFIX]_minInd0.3_fstARRAY.sh
```

#### Weighted FST


``` bash
[PREFIX]_print_idx_minInd0.3_globalFst.sh
```

## Identifying Regions of Divergence associated with Run Timing

Multi-species whole genome Fst scan plot with *esrb* and *lrrc9* highlighted. These plots are the initial run timing phenotype comparisons.

1.  Pink (Even): Early v. Late

2.  Sockeye: Early Creek v. Late Beach

3.  Chum: Early (Summer) v. Late (Fall)

4.  Coho: Early v. Late


``` r
Figure1_Fst_WholeGenomeManhattan_FourSpp.R
```

The two other comparisons were included as a supplementary figure.

5.  Pink (Odd): Early v. Late

6.  Sockeye: Early Creek v. Late Creek


``` r
FigureS2_Fst_WholeGenomeManhattan_OddPink-Sockeye.R
```

The two highlighted regions of divergence that were shared across 2+ species were *lrrc9* and *esrb*. We calculated the mean weighted Fst across the entire genome for each species. We also calculated local weighted Fst within the shared elevated peaks (that have different bounds for each species).


``` r
TableS2_fourspp_weightedFst_calculations.R
```

*Lrrc9* and *esrb* are in the same \~850 kb duplicated region across chromosomes 29 and 35. We plotted the genes in this region to better visualize the duplication and the surrounding regions. To do so, we downloaded the gff file from the reference genome on NCBI. For the zoom in plot, we made the values on chr35 negative to avoid plotting as an inversion.


``` r
Figure2_A_WGD_Illustration.R
```

The gene/duplicate figures were imported into Inkscape and further developed into an illustration for the manuscript.

### Consensus Tree

First, NCBI sequences were taken for *lrrc9*, *esrb*, and their duplicated variants (*lrrc9*-like, and *esrb*-like) from reference genomes in the four salmon species included in this study. Pike was used as an outgroup. They were imported into Geneious Prime and phylogenetic trees were created and exported as newick files.

@PatrickBarry-NOAA - any additional comments?

Consensus trees were plotted for *lrrc9* and *esrb* using the code below.


``` r
Figure2_BC_ConsensusTrees.R
```

### Genotype Assignment from Local PCAs & Resulting Fsts

The peaks varied across species in boundary width. To assign genotypes across species, we used the overlapping peak region as the peak boundary. For *lrrc9,* the shared peak was directly in the *lrrc9* gene region. For *esrb*, it was directly around *esrb* where coho and chum share a peak.

***lrrc9***: NC_068455.1:28128954-28169980

***esrb***: NC_068449.1:25414060-25501622


``` bash
sbatch lrrc9_pca_minInd0.3.sh
sbatch esrb_pca_minInd0.3.sh
```

Plotted each species-gene local PCA independently. In each of the R scripts, a local PCA was plotted, where the PC1 showed variance across run timing phenotypes. Each local PCA was split into three clusters along the PC1 axis. The three groups were defined as the homozygous early (EE), heterozygous (EL), and homozygous late (LL).


``` r
../species_specific/[PREFIX]_lrrc9_pca.R
../species_specific/[PREFIX]_esrb_pca.R
```

The scripts above also split the individuals by sampleID into these three genotypes. Heterozygous individuals were dropped, and EE and LL individuals were placed into separate bamslists. The bamslists were used to run an allele-based FST, which was run in the scripts below for each species comparison.


``` bash
sbatch [PREFIX]_minInd0.3_lrrc9_allele_fst.sh
sbatch [PREFIX]_minInd0.3_esrb_allele_fst.sh
```

These figures are the local PCAs for each species with the cutoff used to delineate the genotypes. Also, the zoomed in FST on the main regions of differentiation using the pca-assigned genotypes. Genes within the regions are also included.


``` r
Figure3_Allele_LocalPCAandFst.R
```

Supplemental local PCAs - adding Whitefish to the local *lrrc9* sockeye PCA, and Pink Odd lineage local PCA. Whitefish individuals were included in the EE and LL bamslists for the above FST.


``` r
FigureS3_AlleleLocalPCAs_OddPink-Sockeye.R
```

For the *esrb* and *lrrc9* regions plotted in Figure 3B&C, we determined if there were SNPs with elevated FST (FST \> 0.5) shared across species. Furthermore, if SNPs were within an annotated genes region, the genes and associated GO terms were listed with that SNP.


``` r
TableS4_allele-based_FST_positions_genes.R
```

#### Barplots

Barplots of putative allele proportions from PCA assignment of genotypes.

For pink-odd, we used the GTseq proportions from 2019 instead of whole genome results.


``` r
FigureS4_Lrrc9_Esrb_Barplots.R
```

### Homozygous Allele Phylogenetic Trees

Use the IBS matrix in angsd to create phylogenetic trees from our GL data. This was used within the *lrrc9* gene region and the *esrb* gene region. The associated new and older flags include:

-   -doCov 1, -makeMatrix 1, & -doIBS 1

-   *minMapQ* [20]: the minimum map quality is 20.

-   *minQ* [20]: the minimum base quality is 20.

-   *doMajorMinor* [4]: Set the major and minor as the reference genome major/minor, was changed from 1.

-   *C* [50]: use an adjustment of the Map Quality for excessive mismatches of 50. This is the suggested value for BWA and leads to fewer false positive variant calls

-   *SNP_pval* [1e-10]: only work with sites with a p-value from the likelihood test of less than 1e-10.

-   *GL* [1]: use SAMtools genotype likelihood calling algorithm. One assumption of the SAMtools algorithm is that it assumes that errors are not independent, but that once a first error occurs at a certain site in an individual, a second error is more likely to occur at the same site

-   *minMaf* [0.05]: only work with sites with a maf above 0.05

### Allele-Based Neighbor Joining Tree

All of the bam files were already created at this point but are stored in the species-specific run timing folders. Therefore, the species bamslists have to be concatenated.


``` bash
cd runtiming/

cat ./pink/pink-chum_filtered_bamslist.txt ./sockeye/sock-chum_filtered_bamslist.txt ./anvil/euclide_filtered_bamslist.txt ./coho/coho-chum_filtered_bamslist.txt ./chum/chumrun_bamslist.txt > ./fourspecies/fourspp_bamslist_all.txt
```

There are a total of 432 individuals.

#### Allele Trees

We used a subset of the total individuals assigned to homozygous early (EE) and late genotypes (LL) to plot a digestible phylogenetic tree for the manuscript (ten per species, 5 per genotype where applicable). Individuals were selected largely based on depth, as determined via the following script:


``` r
Individual_Subset_for_Figure4_AlleleNJTree.R
```

The script results in the following two files:

```         
fourspp_lrrc9_top5_ibs_input.txt
fourspp_esrb_top5_ibs_input.txt
```

Upload the files to Slurm Manager and use as input bamslists in the following shells.


``` bash
top5_lrrc9_fourspp_ibs.sh
top5_esrb_fourspp_ibs.sh
```

Plot Trees in R for both *lrrc9* and *esrb*.


``` r
Figure4_AlleleNJTree.R
```
