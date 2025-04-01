# calculate weighted FST WHOLE GENOME pink, sockeye, chum, coho aligned to chum 
# 03/21/2025

packages_needed <- c("ggplot2", "scales", "ggpubr", "ggrepel", "tidyverse",
                     "stringr", "data.table", "plyr","gtools","reshape2", 
                     "here", "magrittr", "patchwork", "cowplot")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

### Read in each chr IDX FST file #########################################

########## PINK-EVEN
pink_list <- as.list(Sys.glob("./results/fst/*pink-chum*early-late*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
pink_df <- pink_list %>%
  set_names(nm = pink_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Pink-Even")

########## PINK-ODD
odd_list <- as.list(Sys.glob("./results/fst/*pink-odd*early-late*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
odd_df <- odd_list %>%
  set_names(nm = odd_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Pink-Odd")

########## SOCKEYE (Creek)
sock_list <- as.list(Sys.glob("./results/fst/*sock-chum*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
sock_df <- sock_list %>%
  set_names(nm = sock_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Sockeye")

########## EUCLIDE (Early Creek v Late Beach)
eucl_list <- as.list(Sys.glob("./results/fst/*euclide*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
eucl_df <- eucl_list %>%
  set_names(nm = eucl_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Euclide")

########## CHUM
chum_list <- as.list(Sys.glob("./results/fst/*chumrun*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
chum_df <- chum_list %>%
  set_names(nm = chum_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Chum")


########## COHO
coho_list <- as.list(Sys.glob("./results/fst/*coho-chum*minInd0.3*.idx.txt"))

# read in all data files that match wildcard
coho_df <- coho_list %>%
  set_names(nm = coho_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = F, sep = "\t", row.names = NULL, 
                  col.names = c("chrName", "midPos", "A", "B"))
  ) %>% 
  mutate(Species = "Coho")

rm(pink_list, odd_list, sock_list, eucl_list, chum_list, coho_list)

##### Combine each IDX Fst Comparison #######################
# bind rows together
idx_df <- bind_rows(pink_df, odd_df, sock_df, eucl_df, chum_df, coho_df) %>%
  mutate(A = as.numeric(A),
         B = as.numeric(B))

idx_df$A[idx_df$A < 0] <- 0

weighted_Fsts <- idx_df %>%
  group_by(Species) %>%
  dplyr::summarise(globalFst = sum(A)/sum(B))

write.csv(weighted_Fsts, "./results/fst/weighted_globalFSTs.csv",
            row.names = F, quote = F)

rm(pink_df, odd_df, sock_df, eucl_df, chum_df, coho_df) # remove to save space


### Calculate Peak Boundaries with Fst ##########################################

# lrrc9
pink_Fst <- read.delim2("./results/fst/pink-chum_NC_068455.1_early-late_minInd0.3.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Pink-Even")
odd_Fst <- read.delim2("./results/fst/pink-odd_NC_068455.1_early-late_minInd0.3.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Pink-Odd")
eucl_Fst <- read.delim("./results/fst/euclide_NC_068455.1_beach-creek_minInd0.3.sfs.pbs.fst.txt",
                       row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Euclide")
sock_Fst <- read.delim("./results/fst/sock-chum_NC_068455.1_early-late_minInd0.3.sfs.pbs.fst.txt",
                       row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Sockeye")
chumL_Fst <- read.delim2("./results/fst/chumrun_NC_068455.1_fall-summer_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Chum")

#esrb
chumE_Fst <- read.delim2("./results/fst/chumrun_NC_068449.1_fall-summer_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                         row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Chum")
coho_Fst <- read.delim2("./results/fst/coho-chum_NC_068449.1_early-late_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                         row.names = NULL,sep = "\t")[,c(2:3,5)] %>% mutate(Species = "Coho")

# Bind fst dataframes and make all columns as.character
dfs <- list(pink_Fst, odd_Fst, sock_Fst, eucl_Fst, chumL_Fst, chumE_Fst, coho_Fst)

dfs <- lapply(dfs, function(df) {
  df[] <- lapply(df, as.character)  # Convert all columns to character
  return(df)
})

all_fst <- do.call(rbind, dfs)
colnames(all_fst) <- c("chrName", "midPos","Fst", "Species")

rm(pink_Fst, odd_Fst, sock_Fst, eucl_Fst, chumL_Fst, chumE_Fst, coho_Fst, dfs)

# First remove all sites that were less than 0.5
gene_fst <- all_fst %>%
  mutate(midPos = as.numeric(midPos), 
         Fst = as.numeric(Fst)) %>%
  filter(Fst > 0.5)

# define the broad boundaries for chr29 & chr35
gene_fst <- rbind(gene_fst %>%
                    filter(chrName == "NC_068449.1",
                           midPos > 24*1e6, midPos < 27*1e6),
                  gene_fst %>%
                    filter(chrName == "NC_068455.1",
                           midPos > 26*1e6, midPos < 29*1e6))

# Calculate gene boundaries
# Note that if there are no Fst values above 0.5, it will be dropped
gene_boundary <- gene_fst %>%
  group_by(Species, chrName) %>%
  dplyr::summarise(minPos = min(midPos), maxPos = max(midPos)) %>%
  mutate(Boundary = paste0(minPos,"-",maxPos),
         BoundarySize = maxPos - minPos) %>%
  ungroup()
View(gene_boundary)

# Use gene boundaries to define min and max position
# Not that the mutate function was used for species that weren't assigned appropriate boundaries
idx_calcs <- idx_df %>%
  filter(chrName %in% unique(gene_boundary$chrName)) %>%
  dplyr::left_join(gene_boundary, by = c('Species','chrName')) %>%
  mutate(maxPos = case_when(chrName == "NC_068449.1" ~ maxPos,
                            Species == "Pink-Odd" ~ gene_boundary$maxPos[gene_boundary$Species=="Pink-Even"],
                            Species == "Sockeye" ~ gene_boundary$maxPos[gene_boundary$Species=="Euclide"],
                            T ~ maxPos),
         minPos = case_when(chrName == "NC_068449.1" ~ minPos,
                            Species == "Pink-Odd" ~ gene_boundary$minPos[gene_boundary$Species=="Pink-Even"],
                            Species == "Sockeye" ~ gene_boundary$minPos[gene_boundary$Species=="Euclide"],
                            T ~ minPos)) %>%
  filter(midPos >= minPos, midPos <= maxPos)

idx_local_summary <- idx_calcs %>%
  group_by(Species, chrName) %>%
  dplyr::summarize(peak_mean_Fst = sum(A)/sum(B)) %>%
  ungroup() %>%
  left_join(gene_boundary, by = c("Species","chrName")) %>%
  left_join(weighted_Fsts, by = "Species") %>%
  select(Species, globalFst, everything())

write.csv(idx_local_summary, "./results/fst/weighted_global-peak_meanFsts.csv",
          row.names = F, quote = F)
