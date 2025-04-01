# quick plots of lrrc9 only

library(here)
library(tidyverse)

SPECIES = 'Chum'

METADATAFILE <- "./data/raw/fourspecies_runtiming_metadata.csv"
BAMFILE <- "./data/bams/chumrun_bamslist.txt"

########### ALL POPS #################################

# read in the covariance matrix
cov <- as.matrix(read.table(paste0("./results/pca/chumrun_NC_068455.1_lrrc9_minInd0.3_minDepthHalf.cov")))
  e <- eigen(cov)                            # calculate eigenvector values
  e_vectors <- as.data.frame(e$vectors) %>%
    mutate(FID = row_number())               # add FID to e_vectors
  e_per <- e$values/sum(e$values)            # percent explained by each component

# call in bams
bam_df <- read.table(BAMFILE, header = F) %>%
  mutate(ABLG = as.numeric(gsub("[^0-9-]", "", V1))) %>%
  select(ABLG)
head(bam_df)

# call in some metadata
pop_df <- read.csv(METADATAFILE, header = T) %>%
  filter(Species == SPECIES) %>%
  mutate(ABLG = as.numeric(sub('ABLG','',sampleID))) 

popFID <- inner_join(bam_df, pop_df, by = "ABLG")

##combine row names (population info) with the covariance matrix
pca.vectors = as_tibble(cbind(popFID, e_vectors))

# determine the variance explained as a percent
pca.eigenval.sum = sum(e$values) #sum of eigenvalues
  varPC1 <- (e$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
  varPC2 <- (e$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2

#### PLOTTING

# color palette
mypalette <- c("goldenrod1", "royalblue3")
pca.vectors$Runtime <- factor(pca.vectors$Runtime, levels = unique(pca.vectors$Runtime))
names(mypalette) <- levels(pca.vectors$Runtime)

theme_set(
  theme(legend.title=element_blank(), legend.text=element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(angle = 0, size = 10), axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white"), panel.spacing = unit(0,"lines"),
        legend.position = "right", strip.text.y = element_text(angle = 0)
  )
)

ggplot(data = pca.vectors, aes(x=V1, y=V2, fill = Runtime, color = Runtime)) + 
  geom_point(alpha = 0.7, size = 1.8, pch = 21) +
  scale_fill_manual(values = mypalette) +
  scale_color_manual(values = mypalette) +
  ggtitle("Chum: PCA of lrrc9 gene") + 
  labs(x = paste0("PC1 (",round(varPC1, digits = 2),"%)"), 
       y= paste0("PC2 (",round(varPC2, digits = 2),"%)")) 


### FIND INDIVIDUALS ASSOCIATED WITH EACH GENOTYPE ###########################

EE_cutoff <- 0.1
LL_cutoff <- -0.05

ggplot(data = pca.vectors, aes(x=V1, y=V2, fill = Runtime, color = Runtime)) + 
  geom_point(alpha = 0.7, size = 1.8, pch = 21) +
  theme(legend.position = "right") +
  scale_fill_manual(values = mypalette) +
  scale_color_manual(values = mypalette) +
  scale_x_continuous(breaks = c(-0.05, 0, 0.05, 0.1, 0.15)) +
  geom_vline(xintercept = EE_cutoff, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = LL_cutoff, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle("Chum: PCA of lrrc9 gene") + 
  labs(x = paste0("PC1 (",round(varPC1, digits = 2),"%)"), 
       y= paste0("PC2 (",round(varPC2, digits = 2),"%)"))

print(paste0(nrow(pca.vectors %>% filter(V1 < EE_cutoff, V1 > LL_cutoff))," heterozygote individuals"))

allele_groups <- pca.vectors[,1:10] %>% 
  mutate(Genotype = case_when(V1 > EE_cutoff ~ "EE",  
                             V1 < LL_cutoff ~ "LL",
                             TRUE ~ "EL")) %>% 
  mutate(name = str_c("ABLG",ABLG))

# write out to text file to be transferred and used as input
write.table(allele_groups[,c('name','Genotype')], file = paste0("./data/R/chumrun_lrrc9_alleles.txt"), 
            col.names = F, row.names = F, sep = "\t", quote = F)

allele_props <- allele_groups %>%
  group_by(Genotype, Runtime) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(prop = n/nrow(pca.vectors))
allele_props

write.table(allele_props, 
            file = paste0("./data/R/chumrun_lrrc9_allele_proportions.txt"), 
            col.names = T, row.names = F, sep = "\t", quote = F)

