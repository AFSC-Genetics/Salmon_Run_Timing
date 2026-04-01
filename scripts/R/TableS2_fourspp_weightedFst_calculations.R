# calculate weighted FST WHOLE GENOME pink, sockeye, chum, coho aligned to chum 
# 03/21/2025
# Updated for esr1 on 03/19/2026

# setwd("../../Salmon_Run_Timing")

packages_needed <- c("ggplot2", "scales", "ggpubr", "ggrepel", "tidyverse",
                     "stringr", "data.table", "plyr","gtools","reshape2", 
                     "here", "magrittr", "patchwork", "cowplot")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

### Read in each chr IDX FST file #########################################

########## PINK-EVEN
pink_list <- as.list(Sys.glob("./results/fst/pink/idx/pink-chum*early-late*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
pink_df <- pink_list %>%
  set_names(nm = pink_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Pink-Even")

########## PINK-ODD
odd_list <- as.list(Sys.glob("./results/fst/pink/idx/pink-odd*early-late*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
odd_df <- odd_list %>%
  set_names(nm = odd_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Pink-Odd")

########## SOCKEYE (Creek)
sock_list <- as.list(Sys.glob("./results/fst/sockeye/idx/sock-chum*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
sock_df <- sock_list %>%
  set_names(nm = sock_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Sockeye")

########## EUCLIDE (Early Creek v Late Beach)
eucl_list <- as.list(Sys.glob("./results/fst/sockeye/idx/euclide*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
eucl_df <- eucl_list %>%
  set_names(nm = eucl_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Sockeye - Euclide")

########## EUCLIDE (Early Creek v Late Beach)
pick_list <- as.list(Sys.glob("./results/fst/sockeye/idx/pick-chum*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
pick_df <- pick_list %>%
  set_names(nm = pick_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Sockeye - Pick")

########## CHUM
chum_list <- as.list(Sys.glob("./results/fst/chum/idx/chumrun*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
chum_df <- chum_list %>%
  set_names(nm = chum_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Chum")

########## COHO
coho_list <- as.list(Sys.glob("./results/fst/coho/idx/coho-chum*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
coho_df <- coho_list %>%
  set_names(nm = coho_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Coho")

rm(pink_list, odd_list, sock_list, eucl_list, pick_list, chum_list, coho_list)

##### Combine each IDX Fst Comparison #######################
# bind rows together
idx_df <- bind_rows(pink_df, odd_df, 
                    sock_df, eucl_df, pick_df, 
                    chum_df, 
                    coho_df) %>%
  mutate(A = as.numeric(A),
         B = as.numeric(B))

idx_df$A[idx_df$A < 0] <- 0

weighted_Fsts <- idx_df %>%
  group_by(Species) %>%
  dplyr::summarise(mean_Fst = sum(A)/sum(B)) %>%
  arrange(mean_Fst)

write.csv(weighted_Fsts, "./results/fst/weighted_meanFSTs_2026.csv",
            row.names = F, quote = F)

rm(pink_df, odd_df, sock_df, eucl_df, chum_df, coho_df) # remove to save space


### Calculate Peak Boundaries with Fst ##########################################

# lrrc9
pink_Fst <- read.delim2("./results/fst/pink/pink-chum_NC_068455.1_early-late_minInd0.3.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Pink-Even")
odd_Fst <- read.delim2("./results/fst/pink/pink-odd_NC_068455.1_early-late_minInd0.3.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Pink-Odd")
eucl_Fst <- read.delim("./results/fst/sockeye/euclide_NC_068455.1_beach-creek_minInd0.3.sfs.pbs.fst.txt",
                       row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Sockeye - Euclide")
sock_Fst <- read.delim("./results/fst/sockeye/sock-chum_NC_068455.1_early-late_minInd0.3.sfs.pbs.fst.txt",
                       row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Sockeye")
chumL_Fst <- read.delim2("./results/fst/chum/chumrun_NC_068455.1_fall-summer_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Chum")

# esrb
chumE_Fst <- read.delim2("./results/fst/chum/chumrun_NC_068449.1_fall-summer_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                         row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Chum")
coho_Fst <- read.delim2("./results/fst/coho/coho-chum_NC_068449.1_early-late_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                         row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Coho")

# esr1
euclE1_Fst <- read.delim("./results/fst/sockeye/euclide_NC_068428.1_beach-creek_minInd0.3.sfs.pbs.fst.txt",
                       row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Sockeye - Euclide")
pickE1_Fst <- read.delim("./results/fst/sockeye/pick-chum_NC_068428.1_early-late_minInd0.3.sfs.pbs.fst.txt",
                         row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Sockeye - Pick")
cohoE1_Fst <- read.delim2("./results/fst/coho/coho-chum_NC_068428.1_early-late_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Coho")


# Bind fst dataframes and make all columns as.character
dfs <- list(pink_Fst, odd_Fst, sock_Fst, eucl_Fst, chumL_Fst, chumE_Fst, coho_Fst, pickE1_Fst, euclE1_Fst, cohoE1_Fst)

dfs <- lapply(dfs, function(df) {
  df[] <- lapply(df, as.character)  # Convert all columns to character
  return(df)
})

all_fst <- do.call(rbind, dfs)
colnames(all_fst) <- c("chrName", "midPos","Fst", "Species")

rm(pink_Fst, odd_Fst, sock_Fst, eucl_Fst, chumL_Fst, chumE_Fst, coho_Fst, pickE1_Fst, euclE1_Fst, cohoE1_Fst, dfs)

# First remove all sites that were less than 0.5
gene_fst <- all_fst %>%
  mutate(midPos = as.numeric(midPos), 
         Fst = as.numeric(Fst)) %>%
  filter(Fst > 0.5)

# define the broad boundaries for chr8, chr29 & chr35
gene_fst <- rbind(gene_fst %>%
                    filter(chrName == "NC_068428.1",
                           midPos > 23*1e6, midPos < 25*1e6), # get into the broad region of interest,
                  gene_fst %>%
                    filter(chrName == "NC_068449.1",
                           midPos > 24*1e6, midPos < 27*1e6), 
                  gene_fst %>%
                    filter(chrName == "NC_068455.1",
                           midPos > 26*1e6, midPos < 29*1e6))

# remove one rogue SNP in chum (from visual inspection of broad general peak region)
gene_fst <- gene_fst %>%
  filter(!(Species == 'Chum' & midPos == 28357254 & chrName == "NC_068455.1"))

# Calculate gene boundaries
# Note that if there are no Fst values above 0.5, it will be dropped
gene_boundary <- gene_fst %>%
  group_by(Species, chrName) %>%
  dplyr::summarise(minPos = min(midPos), maxPos = max(midPos)) %>%
  mutate(Boundary = paste0(minPos,"-",maxPos),
         BoundarySize = maxPos - minPos) %>%
  ungroup() %>%
  arrange(chrName)
# View(gene_boundary)

gene_boundary_by_spp <- gene_boundary %>%
  dplyr::mutate(spp = sub(" - .*", "", Species)) %>% 
  dplyr::group_by(spp, chrName) %>%
  dplyr::mutate(
    minPos = min(minPos, na.rm = TRUE),
    maxPos = max(maxPos, na.rm = TRUE),
    BoundarySize = maxPos - minPos) %>%
  dplyr::ungroup() %>%
  dplyr::select(Species,chrName,minPos,maxPos)

# Use gene boundaries to define min and max position
# Not that the mutate function was used for species that weren't assigned appropriate boundaries
idx_calcs <- idx_df %>%
  filter(chrName %in% unique(gene_boundary$chrName)) %>%
  dplyr::left_join(gene_boundary_by_spp, by = c('Species','chrName')) %>%
  filter(midPos >= minPos, midPos <= maxPos)

idx_local_summary <- idx_calcs %>%
  dplyr::group_by(Species, chrName) %>%
  dplyr::summarize(peak_mean_Fst = sum(A)/sum(B)) %>%
  dplyr::ungroup() %>%
  left_join(gene_boundary_by_spp, by = c("Species","chrName")) %>%
  left_join(weighted_Fsts, by = "Species") %>%
  filter(!is.na(minPos)) %>%
  dplyr::arrange(chrName) %>%
  dplyr::rename(wg_mean_Fst = mean_Fst) %>%
  dplyr::select(Species, chrName, everything())

idx_local_summary

write.csv(idx_local_summary, "./results/fst/peak_weighted_meanFsts.csv",
          row.names = F, quote = F)

#### Determine the broadest peak range across species & associated size ########

idx_local_summary %>%
  dplyr::group_by(chrName) %>%
  dplyr::summarize(min = min(minPos, na.rm = T), 
                   max = max(maxPos, na.rm = T)) %>%
  dplyr::ungroup() %>%
  mutate(PeakSize = max - min)
