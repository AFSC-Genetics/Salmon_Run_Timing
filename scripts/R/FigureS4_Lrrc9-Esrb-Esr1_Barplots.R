# BARPLOTS FOR LRRC9 AND ESRB ALLELE FREQUENCIES
# the proportions are pulled from the outputs of the species-specific pca_lrrc9 and pca_esrb scripts

packages_needed <- c("ggplot2", "scales", "ggpubr", "tidyverse", "grid", "gridExtra",
                     "lattice", "patchwork", "here", "cowplot", "magrittr", "ggh4x")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

# setwd("../../Salmon_Run_Timing")
outdir="./figures/barplot/"
# outdir="../Salmon_runtiming/2024_fourspecies/figures/barplot/"

### Input Proportions##################

#### LRRC9 ##################

pinkeven_lrrc9 <- read.table(file = "data/R/threespp_lrrc9_alleles/pink-chum_minInd0.3_lrrc9_allele_proportions.txt",
                         header = T, sep = "\t")%>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Pink Even")

pinkodd_lrrc9 <- read.csv("./data/R/threespp_lrrc9_alleles/PinkOdd_2019_GTseq.csv", row.names = NULL) %>%
  mutate(Species = "Pink Odd",
         Runtime = Pheno,
         n = nPinks,
         Genotype = case_when(value == 4004 ~ "EE",
                              value == 2004 ~ "EL",
                              value == 2002 ~ "LL")) %>%
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

#### ESRB ####################
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

#### ER1 ####################
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

# add levels for ordering plots
genes_df$Species <- factor(genes_df$Species, levels = c("Pink Even", "Pink Odd", "Sockeye", "Chum", "Coho"))
genes_df$Runtime <- factor(genes_df$Runtime, levels = c("Early","Late","Early Stream","Late Stream", "Late Beach", 
                                                        "Early Stream (Teal)", "Late Stream (Whitefish)", 
                                                        "Early Stream (Pick)", "Late Stream (Pick)", "Late Beach (Anvil)"))
genes_df$Gene <- factor(genes_df$Gene, levels = c("LRRC9", "ESRB", "ESR1"))
levels(genes_df$Species)
levels(genes_df$Gene)

write.csv(genes_df, "./data/R/AlleleGroup_CountPerSpecies-Genotype_2.csv", 
            row.names = F)

#### LRRC9 Plot  ###############

lrrc9_barplot <- ggplot(filter(genes_df, Gene == "LRRC9"), 
                        aes(x = Runtime, y = n, fill = Genotype)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_nested(~ Species, 
               scales = "free", space = "free",
               labeller = label_wrap_gen(10)) +
  scale_fill_manual(name = "Genotype", values = c("goldenrod2", "mediumseagreen", "royalblue3")) +
  ylab("Proportion") + xlab("Run Timing Phenotype") +
  ggtitle(expression(italic(lrrc9))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 9)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 22),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12, vjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, angle = 0, color = "black"),
        axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20, angle = 90, margin = margin(0,5,0,1)),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(0,"mm"),
        strip.background = element_rect(fill = "gray95", color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        legend.background = element_rect(color = "white")
  )
lrrc9_barplot

#### ESRB Plot  ###############
esrb_barplot <- ggplot(filter(genes_df, Gene == "ESRB"), 
                        aes(x = Runtime, y = n, fill = Genotype)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_nested(~ Species, 
               scales = "free", space = "free",
               labeller = label_wrap_gen(10)) +
  scale_fill_manual(name = "Genotype", values = c("goldenrod2", "mediumseagreen", "royalblue3")) +
  ylab("Proportion") + xlab("Run Timing Phenotype") +
  ggtitle(expression(italic(esrb))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 9)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 22),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12, vjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, angle = 0, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(0,"mm"),
        strip.background = element_rect(fill = "gray95", color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        legend.background = element_rect(color = "white")
  )
esrb_barplot

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
  scale_x_discrete(labels = function(x) str_wrap(x, width = 9)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_text(size = 20), legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 22),
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 12, vjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, angle = 0, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(0,"mm"),
        strip.background = element_rect(fill = "gray95", color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        legend.background = element_rect(color = "white")
  )
esr1_barplot

### Combine Plots ##############
final_barplot <- cowplot::plot_grid(lrrc9_barplot, esrb_barplot, esr1_barplot, 
                                    rel_widths = c(10, 4.8, 10),
                                    nrow = 1, align = "h",
                                    labels = c('A','B','C'), 
                                    label_fontfamily = "ArialMT",
                                    label_size = 30, 
                                    label_colour = "black")
# final_barplot

x.grob <- grid::textGrob("Run Timing Phenotype", 
                   gp=gpar(col="black", fontsize=18))



# PDF and JPG outputs
pdf(file = paste0(outdir,"three_genes_barplot_",format(Sys.Date(),"%Y%m%d"),".pdf"), 
    width = 21, height = 6)
grid.arrange(gridExtra::arrangeGrob(final_barplot, bottom = x.grob))
dev.off()

jpeg(file = paste0("./figures/barplot/fourspp_genes_barplot_",format(Sys.Date(),"%Y%m%d"),".jpg"), 
     width = 33, height = 15, res = 200, units = 'cm')
grid.arrange(gridExtra::arrangeGrob(final_barplot, bottom = x.grob))
dev.off()

