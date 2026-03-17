# quick plots of esr1 only

library(here)
library(stringr)
library(tools)
library(tidyverse)

SPECIES='Sockeye'
prefix <- "pick-chum"
name = "NC_068428.1_s23264160_e23645630_minInd0.3_esr1"
BAMFILE <- "./data/bams/pick-chum_filtered_bamslist.txt"
# allmeta_df <- read.csv("./data/raw/fourspecies_runtiming_metadata.csv", header = T)
pick_df <- read.csv("./data/raw/pick_creek_metadata.csv", header = T)

#### Only Pick Creek #### 

# read in the covariance matrix
cov <- as.matrix(read.table(paste0("./results/pca/",prefix,"_",name,".cov")))
  e <- eigen(cov)
  e_vectors <- as.data.frame(e$vectors)
  # determine the variance explained as a percent
  varPC1 <- (e$values[1]/sum(e$values))*100 # PC1 variance
  varPC2 <- (e$values[2]/sum(e$values))*100 # PC2 variance

# call in metadata
bam_df <- read.table(BAMFILE, header = F) %>%
  dplyr::mutate(temp = basename(file_path_sans_ext(V1)),
                temp = gsub("^[^_]*_", "", temp), # remove everything after 1st underscore
                sampleID = str_extract(temp, "[^_]+")) %>%
  select(sampleID)

sock_meta <- pick_df %>%
  mutate(sampleID = as.character(ABLG)) %>%
  inner_join(bam_df, ., by = "sampleID")

##combine row names (population info) with the covariance matrix
pca.vectors = as_tibble(cbind(sock_meta, e_vectors))

###### PLOTTING

pca.vectors$Runtime <- factor(pca.vectors$Runtime, levels = c("Early", "Late"))
  mypalette <- c("goldenrod1", "royalblue3")
  names(mypalette) <- levels(pca.vectors$Runtime)

theme_set(
  theme(legend.position = "right", 
        legend.title=element_blank(), legend.text=element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(angle = 0, size = 10), axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white"), panel.spacing = unit(0,"lines"),
        strip.text.y = element_text(angle = 0)
  )
)

ggplot(data = pca.vectors, aes(x=V1, y=V2, fill = Runtime, color = Runtime)) + 
  geom_point(alpha = 0.7, size = 1.8, pch = 21) +
  theme(legend.position = "right") +
  scale_fill_manual(values = mypalette) +
  scale_color_manual(values = mypalette) +
  scale_x_continuous(breaks = c(-.2,-.1, 0, 0.1, 0.2)) +
  ggtitle("Sockeye - Pick Creek: Chr 8 - esr1 gene") +
  labs(x = paste0("PC1 (",round(varPC1, digits = 2),"%)"), 
       y= paste0("PC2 (",round(varPC2, digits = 2),"%)")) 

### FIND INDIVIDUALS ASSOCIATED WITH EACH "GENOTYPE"

EE_cutoff <- 0.1
LL_cutoff <- -0.1

# plot pca with cutoffs
ggplot(data = pca.vectors, aes(x=V1, y=V2, fill = Runtime, color = Runtime)) + 
  geom_point(alpha = 0.7, size = 1.8, pch = 21) +
  theme(legend.position = "right") +
  scale_fill_manual(values = mypalette) +
  scale_color_manual(values = mypalette) +
  scale_x_continuous(breaks = c(-.2,-.1, 0, 0.1, 0.2)) +
  geom_vline(xintercept = EE_cutoff, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = LL_cutoff, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle("Sockeye (Euclide): PCA of esr1 gene") + 
  labs(x = paste0("PC1 (",round(varPC1, digits = 2),"%)"), 
       y= paste0("PC2 (",round(varPC2, digits = 2),"%)")) +
  theme(legend.position = "bottom")

print(paste0(nrow(pca.vectors %>% filter(V1 > LL_cutoff, V1 < EE_cutoff))," heterozygous individual(s)"))

allele_groups <- pca.vectors[,1:10] %>% 
  mutate(Genotype = case_when(V1 > EE_cutoff ~ "EE",  
                              V1 < LL_cutoff ~ "LL",
                              TRUE ~ "EL")) %>% 
  rename(name = sampleID)

# check that correct number of hets were assigned
allele_groups %>% count(Genotype)

# list of individuals and their genotype groups
write.table(allele_groups[,c('name', 'Genotype')], file = paste0("./data/R/twospp_esr1_alleles/",prefix,"_esr1_alleles.txt"), 
            col.names = F, row.names = F, sep = "\t", quote = F)

# write out to text file to be transferred and used as input for allele fst
allele_groups %>%
  filter(Genotype != "EL") %>%
  mutate(Genotype = paste0(Genotype,"_er1")) %>%
  select(name,Genotype) %>%
  write.table(file = paste0("./data/R/twospp_esr1_alleles/",prefix,"_esr1_homozygotes.txt"), 
            col.names = F, row.names = F, sep = "\t", quote = F)

####### STACKED BAR PLOT ALLELE PROPORTIONS INPUT

allele_props <- allele_groups %>%
  group_by(Genotype, Runtime) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(prop = n/nrow(pca.vectors))
allele_props

write.table(allele_props , 
            file = paste0("./data/R/twospp_esr1_alleles/",prefix,"_esr1_allele_proportions.txt"), 
            col.names = T, row.names = F, sep = "\t", quote = F)
