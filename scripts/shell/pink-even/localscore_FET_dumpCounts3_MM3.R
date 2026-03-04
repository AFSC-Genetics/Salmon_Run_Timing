# for slurm run of FET from counts data
# angsd output for two groups of interest is needed
# doCounts 1 dumpCounts 3 doMaf 1 doMajorMinor 3
# doMajorMinor 3 uses the major/minor of a sites file, which will simplify the code

# had an issue with the libPaths being set in bash script so set here
mamba_lib_path <- "/home/nhowe/.conda/envs/rockfish/lib/R/library"
print(mamba_lib_path)

# Add the mamba library path to the start of the R library path
.libPaths(c(mamba_lib_path, .libPaths()))

args=commandArgs(trailingOnly = TRUE)

homedir <- args[1]
species_prefix <- args[2] # the prefix from the slurm runs
run_suffix <- args[3] # the suffix from the slurm runs

# Steps
# 1. Convert doCounts output to usable count file with Total, Minor Counts and Major Counts
# 2. Fisher's exact test for every locus resulting in p-value.
# 3. Local score uses p-values to determine peaks.

# ONLY STEPS 1-2 HERE. 

# WHERE TO GET DATA AND WHERE TO WRITE IT OUT
indir <- paste0(homedir,"/diversity/")
outdir <- paste0(homedir,"/localscore/")

### Install Packages #########################################

packages_needed <- c("plyr", "stringr", "ggplot2", "tidyverse", "data.table")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

rm(packages_needed,i)
### Import data #########################################

#### Pos Data ###########################
posfiles <- Sys.glob(paste0(indir,species_prefix,'*',run_suffix,".pos.gz"))
posnames <- str_extract(posfiles, paste0("(?<=",species_prefix,"_).*(?=_",run_suffix,")")) # extract the string prior to the suffix and including the chr number
group <- str_split_i(posnames,"_",3) # select the second element

# read in pos files
pos_input <- lapply(posfiles, function(x) {
  read_tsv(file = x,
           col_names = T, show_col_types = F,
           skip_empty_rows = TRUE, name_repair = "minimal")
}) %>%
  set_names(group)

print(paste("Positional input has",length(pos_input)/2,"files per group"))

#### MAF Data ###########################
mafsfiles <- Sys.glob(paste0(indir,species_prefix,'*',run_suffix,".mafs.gz"))

mafs_input <- lapply(mafsfiles, function(x) {
  read_tsv(file = x,
           col_names = T, show_col_types = F,
           skip_empty_rows = TRUE, name_repair = "minimal")
}) %>%
  set_names(group)

print(paste("MAF input has",length(mafs_input)/2,"files per group"))

#### Counts Data ###########################
countfiles <- Sys.glob(paste0(indir,species_prefix,'*',run_suffix,".counts.gz"))

# read in doCounts output for each group
counts_input <- lapply(countfiles, function(x) {
  read_tsv(file = x,
           col_names = T, show_col_types = F,
           skip_empty_rows = TRUE, name_repair = "minimal")
})

print(paste("Counts input has",length(counts_input)/2,"files per group"))

#### Merge Pos and Counts Data #######################

counts_pos_df <- pmap(list(pos_input, mafs_input, counts_input), 
                      ~bind_cols(..1, ..2, ..3)) %>%
  bind_rows(.id = "group") %>%
  rename(maf = knownEM, 
         A = totA, C = totC, G = totG, `T` = totT) %>%
  select(group,chr,pos,major,minor,maf,totDepth,A,C,G,`T`)

# remove
rm(counts_input, pos_input, mafs_input, posnames, posfiles, mafsfiles, countfiles, group)

print(paste("Before removing triallelic sites, df has",nrow(counts_pos_df),"sites"))

# Remove triallelic sites and address monoallelic sites
biallelic_df <- counts_pos_df %>%
  rowwise() %>%
  mutate(MajCount = get(major),
         MinCount = get(minor),
         sumMinMaj = MinCount + MajCount,
         ) %>%
  ungroup() %>%
  filter(totDepth == sumMinMaj) # if totDepth != major + minor counts, triallelic

print(paste("After removing triallelic sites, df has",nrow(biallelic_df),"sites"))

# relist by group
count_list <- biallelic_df %>%
  select(group,chr,pos,major,minor,MinCount,MajCount,totDepth) %>%
  split(., .$group)

### Merge chr_pos and compare groups ######################################
# inner join will remove sites that were only in one group
count_df <- inner_join(count_list[[1]], 
                       count_list[[2]] %>% select(-major,-minor), 
                       by = c("chr", "pos")) %>%
  select(!starts_with("group"))

print(paste("After removing sites unique to a group, df has",nrow(count_df),"sites"))

# checks
  if(nrow(count_df %>% filter(totDepth.x != MajCount.x + MinCount.x)) > 0){
    print(" WARNING! ERROR! ALLELE 1 + ALLELE 2 != totDepth ALLELE COUNT IN POP 1!")
    quit(save = "no", status = 1)
  }else if(nrow(count_df %>% filter(totDepth.y != MajCount.y + MinCount.y)) > 0){
    print(" WARNING! ERROR! ALLELE 1 + ALLELE 2 != TOTAL ALLELE COUNT IN POP 2!")
    quit(save = "no", status = 1)
  }else{
    print("Passes checks")
  }
  
### FISHERS EXACT TEST #######################################################
# ONLY DOING THAT IF FISHERS TEST ISNT INCLUDED
# FISHERS TEST TAKES AWHILE SO DON'T ACCIDENTALLY RUN IT IF YOU ALREADY HAVE THE OUTPUT!

# remove unnecessary columns for fishers exact test, which will then lead to local score
fet_df <- count_df %>%
  select(-c(totDepth.x, totDepth.y,major,minor))


# conduct fishers exact test on each row (locus)
fet_df$pval <- apply(fet_df, 1, 
                     function(x){
                       tbl <- matrix(as.numeric(x[3:6]), ncol=2, byrow=T)
                       fisher.test(tbl, alternative="two.sided")$p.value
                     })

# remove unnecessary columns for local score & write out
fet_df %>%
  select(chr,pos,pval) %>%
  write_tsv(paste0(outdir,species_prefix,"_",run_suffix,"_postFET.tsv"), 
            col_names = T)

print("FET complete")
print(paste0("stored here: ",outdir,species_prefix,"_",run_suffix,"_postFET.tsv"))
