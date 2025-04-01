# BARPLOTS FOR LRRC9 AND ESRB ALLELE FREQUENCIES
# the proportions are pulled from the outputs of the species-specifc pca_lrrc9 and pca_esrb scripts

packages_needed <- c("ggplot2", "scales", "ggpubr", "tidyverse", "grid", "gridExtra",
                     "lattice", "patchwork", "here", "cowplot", "magrittr", "ggh4x")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

############ LRRC9 INPUT ##################

pinkeven_lrrc9 <- read.table(file = "./data/R/pink-chum_allele_proportions.txt",
                         header = T, sep = "\t")%>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Pink Even")

pinkodd_lrrc9 <- read.csv("./data/R/PinkOdd_2019_GTseq.csv", row.names = NULL) %>%
  mutate(Species = "Pink Odd",
         Runtime = Pheno,
         n = nPinks,
         Genotype = case_when(value == 4004 ~ "EE",
                              value == 2004 ~ "EL",
                              value == 2002 ~ "LL")) %>%
  select(Species, Runtime, Genotype, n)

sock_lrrc9 <- read.table(file = "./data/R/sock-all_lrrc9_allele_proportions.txt",
                         header = T, sep = "\t") %>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Sockeye")

chum_lrrc9 <- read.table(file = "./data/R/chumrun_lrrc9_allele_proportions.txt", 
                         header = T, sep = "\t") %>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Chum")

lrrc9_df <- rbind(pinkeven_lrrc9, pinkodd_lrrc9, sock_lrrc9, chum_lrrc9) %>%
  mutate(Gene = "LRRC9")

######## ESRB ####################
chum_esrb <- read.table(file = "./data/R/chumrun_esrb_allele_proportions.txt", 
                        header = T, sep = "\t") %>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Chum")

coho_esrb <- read.table(file = "./data/R/coho-chum_esrb_allele_proportions.txt", 
                        header = T, sep = "\t") %>%
  select(Runtime, Genotype, n) %>%
  mutate(Species = "Coho")

esrb_df <- rbind(chum_esrb, coho_esrb) %>%
  mutate(Gene = "ESRB")

###################################
# combine esrb and lrrc9 data
genes_df <- rbind(lrrc9_df, esrb_df) %>%
  select(Species, Gene, everything())

# add levels for ordering plots
genes_df$Species <- factor(genes_df$Species, levels = c("Pink Even", "Pink Odd", "Sockeye", "Chum", "Coho"))
genes_df$Runtime <- factor(genes_df$Runtime, levels = c("Early", "Late", "Early Stream", "Late Stream", "Late Beach"))
genes_df$Gene <- factor(genes_df$Gene, levels = c("LRRC9", "ESRB"))
levels(genes_df$Species)
levels(genes_df$Gene)

write.csv(genes_df, "./data/R/AlleleGroup_CountPerSpecies-Genotype.csv", 
            row.names = F)


#############################################################################
# SPLIT LRRC9 AND ESRB and then make A and B
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
        axis.text.x = element_text(size = 14, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 14, angle = 0, color = "black"),
        axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20, angle = 90, margin = margin(0,5,0,1)),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(0,"mm"),
        strip.background = element_rect(fill = "gray95", color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        legend.background = element_rect(color = "white")
  )
lrrc9_barplot

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
  theme(legend.position = "right",
        legend.title = element_text(size = 20), legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 22),
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 14, angle = 0, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(0,"mm"),
        strip.background = element_rect(fill = "gray95", color = "black"),
        strip.text.x = element_text(size = 14, color = "black"),
        legend.background = element_rect(color = "white")
  )
esrb_barplot

final_barplot <- cowplot::plot_grid(lrrc9_barplot, esrb_barplot, rel_widths = c(9, 6.25),
                     nrow = 1, align = "h",
                     labels = c('A','B'), label_fontfamily = "ArialMT",
                     label_size = 30, label_colour = "black")

x.grob <- grid::textGrob("Run Timing Phenotype", 
                   gp=gpar(col="black", fontsize=18))
final_barplot


# PDF and JPG outputs
pdf(file = paste0("./figures/barplot/fourspp_genes_barplot_",format(Sys.Date(),"%Y%m%d"),"_FINAL.pdf"), 
    width = 12, height = 6)
grid.arrange(gridExtra::arrangeGrob(final_barplot, bottom = x.grob))
dev.off()

jpeg(file = paste0("./figures/barplot/fourspp_genes_barplot_",format(Sys.Date(),"%Y%m%d"),".jpg"), 
     width = 33, height = 15, res = 200, units = 'cm')
grid.arrange(gridExtra::arrangeGrob(final_barplot, bottom = x.grob))
dev.off()

