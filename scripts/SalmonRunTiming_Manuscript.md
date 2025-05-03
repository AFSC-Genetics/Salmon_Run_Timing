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

Investigating run timing differentiation between early and late returning salmon 
from multiple regions of Alaska. All species were aligned to chum to compare the 
specified regions on the genome that are differentiated across species.

## Sample Collection

Data from Auke Creek pink salmon (even and odd lineage), Wood River sockeye 
salmon, Yukon River chum salmon, and Skilak River coho salmon. Sample collection 
data was compiled from the raw data in each four species manually with three 
columns: sampleID, Runtime, and Species

``` r
# All samples and sequencing information can be found at:
as.vector(Sys.glob("./data/raw/*"))
```

### Site Map

@PatrickBarry-NOAA

``` r
FigureS1_SiteMap.R
```

## Low-coverage whole genome sequencing pipeline
The analyses of low-coverage sequencing data follow a standard workflow 
which are outlined below with all relevant flags. Sequencing projects 
varied by species and so coverage differed by species. Specific differences
are referenced in the supplemental methods.

### Initial Setup
Species were analyzed individually. All scripts to process the 
low-coverage data and referenced in the preceding 
sections are located in the `./scripts/` directory. To facilitate 
location of relevent files, paths shown to 
execute the code assume `scripts` is the current working directory. Note that
all .sh scripts have directory structure not reflected in the 
github repository and require some standardization on an HPC. For the purpose 
of this markdown document not all code for species specific analyses are provided.
Instead, script commands are given with PREFIX as a placeholder for the 
species-specific analysis performed.
To reproduce results, one will need to substitute 'PREFIX' in the code below for 
one of the following options:

1.  **Chum (Yukon):** "chumrun"
2.  **Coho (Skilak Lake):** "coho-chum"
3.  **Sockeye**:
    i.  **Teal (Early Creek) & Whitefish (Late Creek)**: "sock-chum"
    ii. **Teal (Early Creek) & Anvil (Late Beach)**: "euclide"
    iii. **Teal (Early Creek) & Whitefish (Late Creek) & Anvil (Late Beach)**: "sock-all"
4.  **Pink:**
    i.  **Even**: "pink-chum"
    ii. **Odd**: "pink-odd

For instance, to perform the Pink salmon even lineage analysis one would 
remove `PREFIX` and insert `pink-chum`. 

All species were aligned to the chum salmon reference genome which can be 
downloaded with 

``` bash
wget ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/373/465/GCF_023373465.1_Oket_V2/GCF_023373465.1_Oket_V2_genomic.fna.gz
```

After unzipping the genome, chromosomes can be listed with

``` bash
grep '^>NC' ./GCF_023373465.1_Oket_V2_genomic.fna |  awk '{print $1}' | sed 's/>//' > chromosomes_all.txt
```

Chum salmon have 37 chromosomal pairs as well as the mtDNA genome. For the 
purpose of this analysis we removed the mitochondrial genome from the list of 
chromosomes to analyze.

``` r
sed '/NC_017838.1/d' chromosomes_all.txt > chromosomes.txt
```

Chromosome sizes were extracted with

``` bash
 cut -f1,2 GCF_023373465.1_Oket_V2_genomic.fna.fai | grep '^NC' > chrom_sizes.txt
```

Two shell scripts index the reference genome. The HPC used for these analyses
used a SLURM manager and so commands presented here use SLURM specific submission
commands such as `sbatch`. If your HPC uses a different manager, these commands
will need to be changed.

``` bash
sbatch ./shell/GCF_023373465.1_Oket_V2_genomic_bwa-index_script.sh
sbatch ./shell/GCF_023373465.1_Oket_V2_genomic_fai_script.sh
```

### Step 1: Quality Control

This step is a quality control of the sequencing data (raw fastqs). It runs 
FastQC [@Andrews2010] on the raw sequence data. FastQC is a java based application 
that runs a set of quality control measures.

Run two shell scripts PREFIX-raw_fastqcARRAY.sh and 
PREFIX-raw_multiqcSLURM.sh can be run in succession.

``` bash
sbatch ./shell/PREFIX-raw_fastqcARRAY.sh
sbatch ./shell/PREFIX-raw_multiqcSLURM.sh
```

### Step 2: Trim

After initial quality control steps, adapters were trimmed from raw fastqs 
with the program TRIMMOMATIC [@Bolger2014] and 
an additional quality-check of the trimmed fastqs with multiQC [@Ewels2016]
was performed. Three shell scripts: PREFIX_trimARRAY.sh, PREFIX-trim_fastqcARRAY.sh, and 
PREFIX-trim_multiqcSLURM.sh should all run in succession.

``` bash
# trim adapters
sbatch ./shell/PREFIX_trimARRAY.sh

# multiqc report
sbatch ./shell/PREFIX-trim_fastqcARRAY.sh
sbatch ./shell/PREFIX-trim_multiqcSLURM.sh
```

The multiQC.html output file was downloaded to evaluate if any issues were 
present in the trimmed data before proceeding to alignment.

### Step 3: Align

Alignment of trimmed reads to the chum salmon reference genome
was performed with BWA and aligned reads were then processed with 
SAMTOOLS [@Li2009]. The average read depth for each sample was then calculated.

``` bash
# align to genome
sbatch ./shell/PREFIX_alignARRAY.sh

# calculate average depths
sbatch ./shell/PREFIX_depthsARRAY.sh
```

The output file PREFIX_depths.csv was imported and run through 
the downsampling script to see if there was a similar distribution of depths across 
early and late individuals. We generally used 0.5x as a depth cutoff; 
removing individual samples which fell below this threshold; however, two chum 
samples were close to the 0.5x cutoff (~0.4x) which were retained.

``` r
./R/species-specific/PREFIX_depth_and_downsampling.R
```

Individuals that failed to pass this filter threshold can be added to a 
file that list individuals to exclude in further analyses (-b flag in 
subsequent scripts). 

#### Mapped Reads

Calculate the total number of mapped reads from the bamslists with samtools.

``` bash
sbatch ./R/PREFIX_mapped_reads.sh
```

The average and standard deviation for the mapped reads after 
downloading the file can be calculated in R with

``` r
mean(read.csv("./results/depth/PREFIX_mapped_reads.csv", header = F)$V1)
sd(read.csv("./results/depth/PREFIX_mapped_reads.csv", header = F)$V1)
```

### Step 4: Genotype Likelihoods

Genotype likelihoods (GL) for polymorphic sites across the genome were calculated 
with ANGSD. This was done by calculating GLs for each chromosome. This script 
uses PREFIX_angsdARRAY_input.txt (array file with chromosomes) and a list of the 
.bam files to analyze (PREFIX_filtered_bamslist.txt).

``` bash
sbatch ./shell/PREFIX_minInd0.3_glsARRAY.sh 
```

Filters slightly differed depending on the species alignment and average coverage.

-   *minMapQ* [15 or 20]: the minimum map quality was 20 for chum salmon and 
     15 for all other species (pink, sockeye, and coho) to relax alignment 
     stringency for those species to the chum reference genome.
-   *minQ* [20]: the minimum quality of nucleotide call was 20 for all species.
-   *C* [50]: An adjustment of the Map Quality for excessive mismatches was 
    set at 50. This is the suggested value for BWA and leads to fewer false positive 
    variant calls.
-   *minDepth* [N or 0.5\*N]: The minimum depth was set to the number of individuals 
    if average depth \> 1 and to half the number of individuals if avg depth 
    was \< 1. Discard site if total sequencing depth (all individuals added 
    together) is below minDepth.
-   *setMaxDepth* [20\*minDepth]: Maximum depth was set to the number of 
    individuals times 10 or 20. This discards the site if total sequencing depth (all 
    individuals added together) is above the number of individuals multiplied 
    by 10 or 20 depending on the average depth of coverage.
-   *SNP_pval* [1e-10]: Retain sites with a p-value from the likelihood 
    test of less than 1e-10.
-   *Gl* [1]: use SAMtools genotype likelihood calling algorithm.
-   *minMaf* [0.05]: Retain sites with a maf above 0.05.
-   *minInd* [0.3\*N]: Set to 30% of sample size to include a proportion of 
    individuals.
-   Additional Filters: *only_proper_pairs* [1]; *remove_bads* [1]; *uniqueOnly* [1]; 
    *trim* [0]; *doMajorMinor* [1]

### Step 5: Collate

The GL output are split by chromosome, which can be combined 
to make a GL file for the whole genome by concatenating the individual beagle 
and maf files. Additionally, a *sites* file can be created that will be used 
as input for $F_{ST}$ calculations under the -sites flag.

``` bash
sbatch ./shell/PREFIX_minInd0.3_concatenate_mafs.sh
sbatch ./shell/PREFIX_minInd0.3_concatenate_beagles.sh
```

### Step 6: Population Analyses - $F_{ST}$

A file containing a list of .bam files (bamslists) for the early and late run 
timing phenotypes (subset of the filtered bamslist of all samples for a given
species) was created to evaluate divergence between run timing phenotypes.

First, angsd calls GLs for each subgroup, but it uses the *sites* flag from the 
sites file created in Step 5, so the SNPs being called will only be from the 
polymorphic files from Step 4. Therefore, some of the filtering steps have been 
removed, such as: *SNP_pval* and *minMaf*. Filters were adjusted for the sample 
size of the subset group: *setminDepth* and set*maxDepth*, and *minInd.*

``` bash
sbatch PREFIX_minInd0.3_fstARRAY.sh
```

#### Weighted $F_{ST}$

``` bash
PREFIX_print_idx_minInd0.3_globalFst.sh
```

## Identifying Regions of Divergence associated with Run Timing

Multi-species whole genome Fst scan plot with *esrb* and *lrrc9* highlighted. 
These plots are the initial run timing phenotype comparisons.

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

The two highlighted regions of divergence that were shared across 2+ species 
were *lrrc9* and *esrb*. We calculated the mean weighted Fst across the entire 
genome for each species. We also calculated local weighted Fst within the shared 
elevated peaks (that have different bounds for each species).

``` r
TableS2_fourspp_weightedFst_calculations.R
```

*Lrrc9* and *esrb* are in the same \~850 kb duplicated region across chromosomes 
29 and 35. We plotted the genes in this region to better visualize the 
duplication and the surrounding regions. To do so, we downloaded the gff file 
from the reference genome on NCBI. For the zoom in plot, we made the values on 
chr35 negative to avoid plotting as an inversion.

``` r
Figure2_A_WGD_Illustration.R
```

The gene/duplicate figures were imported into Inkscape and further developed 
into an illustration for the manuscript.

### Consensus Tree

First, NCBI sequences were taken for *lrrc9*, *esrb*, and their duplicated 
variants (*lrrc9*-like, and *esrb*-like) from reference genomes in the four 
salmon species included in this study. Pike was used as an outgroup. They were 
imported into Geneious Prime and phylogenetic trees were created and exported 
as newick files.

@PatrickBarry-NOAA - any additional comments?

Consensus trees were plotted for *lrrc9* and *esrb* using the code below.

``` r
Figure2_BC_ConsensusTrees.R
```

### Genotype Assignment from Local PCAs & Resulting Fsts

The peaks varied across species in boundary width. To assign genotypes across 
species, we used the overlapping peak region as the peak boundary. For *lrrc9,* 
the shared peak was directly in the *lrrc9* gene region. For *esrb*, it was 
directly around *esrb* where coho and chum share a peak.

***lrrc9***: NC_068455.1:28128954-28169980

***esrb***: NC_068449.1:25414060-25501622

``` bash
sbatch lrrc9_pca_minInd0.3.sh
sbatch esrb_pca_minInd0.3.sh
```

Plotted each species-gene local PCA independently. In each of the R scripts, 
a local PCA was plotted, where the PC1 showed variance across run timing 
phenotypes. Each local PCA was split into three clusters along the PC1 axis. 
The three groups were defined as the homozygous early (EE), heterozygous (EL), 
and homozygous late (LL).

``` r
../species_specific/PREFIX_lrrc9_pca.R
../species_specific/PREFIX_esrb_pca.R
```

The scripts above also split the individuals by sampleID into these three 
genotypes. Heterozygous individuals were dropped, and EE and LL individuals 
were placed into separate bamslists. The bamslists were used to run an 
allele-based FST, which was run in the scripts below for each species comparison.

``` bash
sbatch PREFIX_minInd0.3_lrrc9_allele_fst.sh
sbatch PREFIX_minInd0.3_esrb_allele_fst.sh
```

These figures are the local PCAs for each species with the cutoff used to 
delineate the genotypes. Also, the zoomed in FST on the main regions of 
differentiation using the pca-assigned genotypes. Genes within the regions 
are also included.

``` r
Figure3_Allele_LocalPCAandFst.R
```

Supplemental local PCAs - adding Whitefish to the local *lrrc9* sockeye PCA, 
and Pink Odd lineage local PCA. Whitefish individuals were included in the EE 
and LL bamslists for the above FST.

``` r
FigureS3_AlleleLocalPCAs_OddPink-Sockeye.R
```

For the *esrb* and *lrrc9* regions plotted in Figure 3B&C, we determined if 
there were SNPs with elevated $F_{ST}$ ($F_{ST}$ \> 0.5) shared across species. 
Furthermore, if SNPs were within an annotated genes region, the genes and 
associated GO terms were listed with that SNP.

``` r
TableS4_allele-based_FST_positions_genes.R
```

#### Barplots

Barplots of putative allele proportions from PCA assignment of genotypes.

For pink-odd, we used the GTseq proportions from 2019 instead of whole genome 
results.

``` r
FigureS4_Lrrc9_Esrb_Barplots.R
```

### Homozygous Allele Phylogenetic Trees

Use the IBS matrix in angsd to create phylogenetic trees from our GL data. This 
was used within the *lrrc9* gene region and the *esrb* gene region. The 
associated new and older flags include:

-   -doCov 1, -makeMatrix 1, & -doIBS 1
-   *minMapQ* [20]: the minimum map quality is 20.
-   *minQ* [20]: the minimum base quality is 20.
-   *doMajorMinor* [4]: Set the major and minor as the reference genome 
    major/minor, was changed from 1.
-   *C* [50]: use an adjustment of the Map Quality for excessive mismatches 
    of 50. This is the suggested value for BWA and leads to fewer false positive 
    variant calls
-   *SNP_pval* [1e-10]: only work with sites with a p-value from the likelihood 
    test of less than 1e-10.
-   *GL* [1]: use SAMtools genotype likelihood calling algorithm. One assumption 
    of the SAMtools algorithm is that it assumes that errors are not independent, 
    but that once a first error occurs at a certain site in an individual, a 
    second error is more likely to occur at the same site
-   *minMaf* [0.05]: only work with sites with a maf above 0.05

### Allele-Based Neighbor Joining Tree

All of the bam files were already created at this point but are stored in the 
species-specific run timing folders. Therefore, the species bamslists have to be 
concatenated.

``` bash
cd runtiming/

cat ./pink/pink-chum_filtered_bamslist.txt ./sockeye/sock-chum_filtered_bamslist.txt ./anvil/euclide_filtered_bamslist.txt ./coho/coho-chum_filtered_bamslist.txt ./chum/chumrun_bamslist.txt > ./fourspecies/fourspp_bamslist_all.txt
```

There are a total of 432 individuals.

#### Allele Trees

We used a subset of the total individuals assigned to homozygous early (EE) and 
late genotypes (LL) to plot a digestible phylogenetic tree for the manuscript 
(ten per species, 5 per genotype where applicable). Individuals were selected 
largely based on depth, as determined via the following script:

``` r
Individual_Subset_for_Figure4_AlleleNJTree.R
```

The script results in the following two files:

```         
fourspp_lrrc9_top5_ibs_input.txt
fourspp_esrb_top5_ibs_input.txt
```

Upload the files to Slurm Manager and use as input bamslists in the following 
shells.

``` bash
top5_lrrc9_fourspp_ibs.sh
top5_esrb_fourspp_ibs.sh
```

Plot Trees in R for both *lrrc9* and *esrb*.

``` r
Figure4_AlleleNJTree.R
```
