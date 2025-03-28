# quick plots of lrrc9 only

library(here)
library(stringr)
library(tools)
library(tidyverse)

prefix <- "pink-chum"
suffix <- "minInd0.3"
METADATAFILE <- "./data/raw/pinkEven_collection_updated_metadata.csv"

########### ALL POPS #################################

# read in the covariance matrix
cov <- as.matrix(read.table(paste0("./results/pca/",prefix,"_NC_068455.1_lrrc9_",suffix,".cov")))

# calculate eigenvector values
e <- eigen(cov)
e_values <- e$values
e_vectors <- as.data.frame(e$vectors)

# add FID to e_vectors
e_vectors <- e_vectors %>%
  mutate(FID = row_number())

# percent explained by each component
e_per <- e$values/sum(e$values)

# match populations to data
# call in bams
bam_df <- read.table(paste0("./sedna_files/bams/",prefix,"_filtered_bamslist.txt"), header = F)
bam_df$V1 <- basename(file_path_sans_ext(bam_df$V1))

# convert bam to FID with ABLG
FID <- bam_df %>%
  dplyr::mutate(FID = row_number(),
                temp = gsub("^[^_]*_", "", V1),             # remove everything after 1st underscore
                sampleID = str_extract(temp, "[^_]+")) %>%
  subset(select = -c(V1, temp))

# call in some metadata
pop_df <- read.csv(METADATAFILE, header = T)

# join those two dataframes
popFID <- inner_join(FID, pop_df, by = "sampleID") %>%
  subset(select = -c(FID)) # for some reason pca.vectors wouldn't work with FID column, so removed

##combine row names (population info) with the covariance matrix
pca.vectors = as_tibble(cbind(popFID, e_vectors))

# determine the variance explained as a percent
pca.eigenval.sum = sum(e$values) #sum of eigenvalues
varPC1 <- (e$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (e$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (e$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (e$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

#################
# PLOTTING

# read in population metadata
Runtime <- c("early", "late")
color <- c("goldenrod2", "blue")
df_popMetaData <- data.frame(Runtime, color)
head(df_popMetaData)

mypalette <- df_popMetaData$color
pca.vectors$Runtime <- factor(pca.vectors$Runtime, levels = unique(df_popMetaData$Runtime))

# assign name of gene to the color so it can be plotted properly
names(mypalette) <- levels(pca.vectors$Runtime)
mypalette

theme_set(
  theme( 
    legend.title=element_blank(), 
    legend.text=element_text(size=10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(angle = 0, size = 10),
    axis.title = element_text(size = 12),
    panel.background = element_rect(fill = "white"), 
    panel.spacing = unit(0,"lines"),
    strip.text.y = element_text(angle = 0)
  )
)

pca1_2 <- ggplot(data = pca.vectors, aes(x=V1, y=V2, fill = Runtime, color = Runtime)) + 
  geom_point(alpha = 0.7, size = 1.8, pch = 21) +
  theme(legend.position = "right") +
  scale_fill_manual(values = mypalette) +
  scale_color_manual(values = mypalette) +
  ggtitle(paste0(prefix," chr 35 - lrrc9 gene")) +
  labs(x = paste0("PC1 (",round(varPC1, digits = 2),"%)"), 
       y= paste0("PC2 (",round(varPC2, digits = 2),"%)")) 
pca1_2

ggsave(filename = paste0("./figures/pca/",prefix,"_",suffix,"_lrrc9_pca1_2.jpeg"), plot = pca1_2,
       width = 8, height = 6)

##### CREATE ALLELE CUTOFFS

LL_cutoff <- 0.07
EE_cutoff <- -0.05

# plot pca with cutoffs
pca_cutoff <- ggplot(data = pca.vectors, aes(x=V1, y=V2, fill = Runtime, color = Runtime)) + 
  geom_point(alpha = 0.7, size = 1.8, pch = 21) +
  theme(legend.position = "right") +
  scale_fill_manual(values = mypalette) +
  scale_color_manual(values = mypalette) +
  scale_x_continuous(breaks = c(-0.05, 0, 0.05, 0.1, 0.15)) +
  geom_vline(xintercept = EE_cutoff, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = LL_cutoff, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle("Pink: PCA of LRRC9 gene") + 
  labs(x = paste0("PC1 (",round(varPC1, digits = 2),"%)"), 
       y= paste0("PC2 (",round(varPC2, digits = 2),"%)")) +
  theme(legend.position = "bottom")
pca_cutoff

ggsave(filename = paste0("./figures/pca/lrrc9/",prefix,"_",suffix,"_lrrc9_pca_allele_cutoff.jpeg"), plot = pca_cutoff,
       width = 5, height = 5)

print(paste0(nrow(pca.vectors %>% filter(V1 < LL_cutoff, V1 > EE_cutoff))," heterozygote individuals"))

allele_groups <- pca.vectors[,1:10] %>% 
  mutate(Genotype = case_when(V1 < EE_cutoff ~ "EE",  
                              V1 > LL_cutoff ~ "LL",
                              TRUE ~ "EL")) %>% 
  mutate(name = str_c("ABLG",ABLG)) %>%
  subset(select = c(name, Genotype))

# write out to text file to be transferred and used as input
write.table(allele_groups, file = paste0("../2024_fourspecies/data/R/threespp_lrrc9_alleles/",prefix,"_",suffix,"_lrrc9_alleles_all.txt"), 
            col.names = F, row.names = F, sep = "\t", quote = F)
  # manually added to threespp allele doc in fourspecies project

# remove hets for fst
hom_indivs <- allele_groups %>%
  filter(Genotype != "EL") 

# write out to text file to be transferred and used as input
write.table(hom_indivs , file = paste0("./sedna_files/",prefix,"_",suffix,"_lrrc9_allele_groups.txt"), 
            col.names = F, row.names = F, sep = "\t", quote = F)

########## CREATE STACKED BARGRAPH FOR EARLY AND LATE ALLELE PROPORTIONS

allele_props <- allele_groups %>%
  mutate(Runtime = ifelse(Runtime == "early", "Early", "Late")) %>%
  group_by(Genotype, Runtime) %>%
  summarize(n = n()) %>%
  mutate(prop = n/nrow(pca.vectors))

pink_barplot <- ggplot(allele_props, aes(x = Runtime, y = n, fill = Genotype)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("goldenrod2", "lightpink3", "brown4")) +
  theme_bw() +
  xlab("Pink") +
  theme(legend.position = "right",
        panel.grid.major.y = element_line(color = "gray95"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10, vjust = 0.5),
        axis.text.y = element_text(angle = 0, size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        #axis.title.y = element_text(size = 13, angle = 90, margin = margin(0,5,0,1)),
        title = element_text(size = 12),
        panel.background = element_rect(fill = "white"), 
        panel.spacing = unit(0,"lines"),
        strip.background = element_rect(fill = "gray95"),
        strip.text.y = element_text(angle = 0, size = 12) 
  )
pink_barplot

jpeg(file = paste0("./figures/pca/",prefix,"_barplot_",suffix,".jpg"), 
     width = 3, height = 4, res = 150, units = "in")
print(pink_barplot)
dev.off()

write.table(allele_props , 
            file = paste0("./data/R/",prefix,"_",suffix,"_lrrc9_allele_proportions.txt"), 
            col.names = T, row.names = F, sep = "\t", quote = F)
