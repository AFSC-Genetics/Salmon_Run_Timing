# BARPLOTS FOR LRRC9, ESR1, AND ESRB ALLELE FREQUENCIES
# the proportions are pulled from the outputs of the species-specific pca_lrrc9 and pca_esrb scripts
# Except for the pink odd proportions, which are from GTseq data

packages_needed <- c("ggplot2", "scales", "ggpubr", "tidyverse", "grid", "gridExtra",
                     "lattice", "patchwork", "here", "cowplot", "magrittr", "ggh4x",
                     "svglite")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

outdir="./figures/barplot/"
dir.create(outdir, showWarnings = FALSE)

### Input Proportions ##################

#### LRRC9 ------------------------------------------------  ##################

pinkeven_lrrc9 <- read.table(file = "data/R/threespp_lrrc9_alleles/pink-chum_minInd0.3_lrrc9_allele_proportions.txt",
                             header = T, sep = "\t")%>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Pink Even")

pinkodd_lrrc9 <- read.csv("./results/GTseq/Pink/PinkSalmon_2019Q1_Q4.csv", row.names = NULL) %>%
  mutate(Species = "Pink Odd",
         Runtime = Pheno,
         n = nInd,
         Genotype = case_when(Geno == 'AA' ~ "EE",
                              Geno == 'AT' ~ "EL",
                              Geno == 'TT' ~ "LL")) %>%
  select(Species, Runtime, Genotype, n)

sock_lrrc9 <- read.table(file = "./data/R/threespp_lrrc9_alleles/sock-all_lrrc9_allele_proportions.txt",
                         header = T, sep = "\t") %>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Sockeye")

chum_lrrc9 <- read.table(file = "./data/R/threespp_lrrc9_alleles/chumrun_minInd0.3_minDepthHalf_lrrc9_allele_proportions.txt", 
                         header = T, sep = "\t") %>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Chum")

lrrc9_df <- rbind(pinkeven_lrrc9, pinkodd_lrrc9, sock_lrrc9, chum_lrrc9) %>%
  mutate(Gene = "LRRC9")

#### ESR2B ----------------------------------------------- ####################
chum_esrb <- read.table(file = "./data/R/twospp_esrb_alleles/chumrun_NC_068449.1_s25414060_e25501622_esrb_minInd0.3_minDepthHalf_allele_proportions.txt", 
                        header = T, sep = "\t") %>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Chum")

coho_esrb <- read.table(file = "./data/R/twospp_esrb_alleles/coho-chum_esrb_allele_proportions.txt", 
                        header = T, sep = "\t") %>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Coho")

esrb_df <- rbind(chum_esrb, coho_esrb) %>%
  mutate(Gene = "ESRB")

#### ESR1 ------------------------------------------------ ####################
sock_esr1 <- read.table(file = "./data/R/twospp_esr1_alleles/pick-all_esr1_allele_proportions.txt", 
                        header = T, sep = "\t") %>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Sockeye",
         Runtime = gsub("- ","Stream (",Runtime) %>% gsub("Pick Creek","Pick)",.),
         Runtime = case_when(Runtime=="Late Beach"~'Late Beach (Anvil)',
                             Runtime=="Late Stream" ~ 'Late Stream (Whitefish)',
                             Runtime=="Early Stream" ~ 'Early Stream (Teal)',
                             T ~ Runtime))

coho_esr1 <- read.table(file = "./data/R/twospp_esr1_alleles/coho-chum_esr1_allele_proportions.txt", 
                        header = T, sep = "\t") %>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Coho")

esr1_df <- rbind(sock_esr1, coho_esr1) %>%
  mutate(Gene = "ESR1")

### Combine Gene-Species Proportions #############
genes_df <- rbind(lrrc9_df, esrb_df, esr1_df) %>%
  select(Species, Gene, everything())

genes_df <- genes_df %>%
  mutate(Runtime = gsub('Stream','Crk.',Runtime))

RuntimeList <- c("Early","Late","Early Crk.","Late Crk.", "Late Beach", 
                 "Early Crk. (Teal)", "Late Crk. (Whitefish)", 
                 "Early Crk. (Pick)", "Late Crk. (Pick)", "Late Beach (Anvil)")

# add levels for ordering plots
genes_df$Species <- factor(genes_df$Species, levels = c("Pink Even", "Pink Odd", "Sockeye", "Chum", "Coho"))
genes_df$Gene <- factor(genes_df$Gene, levels = c("LRRC9", "ESRB", "ESR1"))
levels(genes_df$Species)
levels(genes_df$Gene)

write.csv(genes_df, "./data/R/AlleleGroup_CountPerSpecies-Genotype_2.csv", 
          row.names = F)

# expand sockeye location names
genes_df <- genes_df %>%
  mutate(Runtime = Runtime %>% 
           gsub("Early Crk.$",'Early Crk. (Teal)',.) %>% 
           gsub("Late Crk.$",'Late Crk. (White)',.) %>%
           gsub("Whitefish","White",.) %>%
           gsub("Late Beach$",'Late Beach (Anvil)',.)) 

genes_df$Runtime <- factor(genes_df$Runtime, 
                           levels = RuntimeList %>% gsub("fish","",.)) # shorted whitefish to white
levels(genes_df$Runtime)

genes_df$Runtime <- factor(genes_df$Runtime, 
                           levels = c(unique(genes_df$Runtime)[-5],unique(genes_df$Runtime)[5]))

### Plotting ------------------------------------------------------------ ###

# for science, consistent font sizing
width2col <- 7.25  # 2 col, in inches

largeFS <- 10 # gene labels
mediumFS <- 8  # Axis titles, legend titles
smallFS <- 6  # Axis text, annotations, captions

theme_set(
  theme(plot.title = element_text(size = largeFS, hjust = 0.5, margin=margin(5,0,5,0, unit = 'pt')),
        axis.text.x = element_text(size = smallFS, vjust = 0.95, color = "black", 
                                   margin=margin(1,0,0,0, unit = 'pt')),
        axis.text.y = element_text(size = smallFS, angle = 0, color = "black"),
        axis.title.y = element_text(size = mediumFS, angle = 90, margin = margin(0,5,0,1)),
        strip.text.x = element_text(size = mediumFS, color = "black"),
        legend.background = element_rect(color = "white"),
        legend.title = element_text(size = mediumFS), 
        legend.text = element_text(size = mediumFS),
        plot.margin          = margin(t = 0, r = 2, b = 0, l = 2, unit = "pt"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(0,"pt"),
        strip.background = element_rect(fill = "gray95", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
  )
)

#### LRRC9 Plot ###############

lrrc9_barplot <- ggplot(filter(genes_df, Gene == "LRRC9"), 
                        aes(x = Runtime, y = n, fill = Genotype)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_nested(~ Species, 
               scales = "free", space = "free",
               labeller = label_wrap_gen(10)) +
  scale_fill_manual(name = "Genotype", values = c("goldenrod2", "mediumseagreen", "royalblue3")) +
  xlab("Run Timing Phenotype") +
  ggtitle(expression(italic(lrrc9))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme(legend.position = "right", 
        axis.title.y = element_blank())

# lrrc9_barplot

#### ESRB Plot  ###############
esrb_barplot <- ggplot(filter(genes_df, Gene == "ESRB"), 
                       aes(x = Runtime, y = n, fill = Genotype)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_nested(~ Species, 
               scales = "free", space = "free",
               labeller = label_wrap_gen(10)) +
  scale_fill_manual(name = "Genotype", values = c("goldenrod2", "mediumseagreen", "royalblue3")) +
  xlab("Run Timing Phenotype") +
  ggtitle(expression(italic(esr2b))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

# esrb_barplot

#### ESR1 Plot  ###############
esr1_barplot <- ggplot(filter(genes_df, Gene == "ESR1"), 
                       aes(x = Runtime, y = n, fill = Genotype)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_nested(~ Species, 
               scales = "free", space = "free",
               labeller = label_wrap_gen(10)) +
  scale_fill_manual(name = "Genotype", values = c("goldenrod2", "mediumseagreen", "royalblue3")) +
  ylab("Proportion") + xlab("Run Timing Phenotype") +
  ggtitle(expression(italic(esr1))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme(legend.position = "none",
        axis.title.x = element_blank())

# esr1_barplot

### Combine Plots ##############

final_plot <- (esr1_barplot | esrb_barplot | lrrc9_barplot) +
  plot_layout(widths = c(7, 4.5, 10)) +
  plot_annotation(tag_levels = 'A',
                  caption = "Run Timing Phenotype") & 
  theme(
      plot.tag = element_text(size = largeFS, family = "Arial", face = "bold"),
      plot.tag.position = c(0.04, 0.95),
      plot.caption = element_text(size = mediumFS, family = "Arial", hjust = 0.5,
                                  margin = margin(0,0,2,0)),
      legend.box.spacing = unit(0, "pt"),
      legend.box.margin  = margin(0, 2, 0, 2, unit = "pt"),
      legend.margin = margin(0, 1, 0, 1, unit = "pt")
  )

# for high quality svg image
svglite(str_c(outdir,"three_genes_barplot_patchwork.svg"),
        width=width2col, height=width2col/2.5)
final_plot
dev.off()

# for putting directly into manuscript
jpeg(str_c(outdir,"three_genes_barplot_patchwork.jpeg"),
     width=width2col, height=width2col/2.5, res=300, units = "in")
final_plot
dev.off()

