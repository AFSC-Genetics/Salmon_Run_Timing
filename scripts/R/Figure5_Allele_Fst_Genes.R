# COMBINE ALLELE FSTs for LRRC9, ESRB, and ESR1 
# FIGURE 4:
#  A) LRRC9 ALLELE FST
#  B) ESRB ALLELE FST
#  C) ESR1 ALLELE FST

# for science, consistent font sizing
width2col <- 7.25  # 2 col, in inches

largeFS <- 10 # Panel Labels
mediumFS <- 8  # Axis titles, legend titles
smallFS <- 7  # Axis text, annotations, captions

pointSize = 1

# color blind color palette change
# esr1_col = "#009E73"
# esrb_color = "#56B4E9"
# lrrc9_color = "#CC79A7"

# Note: Annotated (characterized) genes for regions with high FST are plotted. 
#  Uncharacterized genes & genes without high FST SNPs (Fst > 0.5) were
#   grayed out and not labelled in the legend to focus on the genes of interest.
packages_needed <- c("ggplot2", "scales", "ggpubr", "ggrepel", "stringr", 
                     "data.table", "plyr","tools","gtools","reshape2", 
                     "patchwork", "cowplot", "tidyverse", "svglite")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

rm(i, packages_needed)

allmeta_df <- read.csv("./data/raw/fourspecies_runtiming_metadata.csv", header = T)

### GFF NCBI Data ###########################################################

# find exons from gff file for genes of interest (from NCBI chum reference genome: GCF_023373465.1)
# gff_df <- read.delim('../Salmon_runtiming/2024_fourspecies/data/R/genomic.gff', header = F, comment.char = "#")

gff_df <- read.delim('./data/R/genomic.gff', header = F, comment.char = "#")
gff_df <- gff_df[,c(1:5,9)] # remove excess columns
colnames(gff_df) <- c("chrName", "RefSeq","exon","start.pos","fin.pos", "ID")

# only keep chr29 & 35 & 8
gff_chr35 <- gff_df %>%
  filter(chrName == "NC_068455.1")
gff_chr29 <- gff_df %>%
  filter(chrName == "NC_068449.1")
gff_chr8 <- gff_df %>%
  filter(chrName == "NC_068428.1")
rm(gff_df)

# prep pattern for str_match below
gene_pattern <- "gene=\\s*(.*?)\\s*;"           # keep string btwn "gene=" & ":product" 
exon_pattern <- "ID=exon-\\s*(.*?)\\s*;Parent"  # keep string btwn "exon=" & ";Parent" 
descr_pattern <- ";description=\\s*(.*?)\\s*;"  # keep description string

### 2A) LRRC9 FST #############################################################

# Define Boundaries
xstart.lrrc9 = 27.86
xend.lrrc9 = 28.24
pca.start.lrrc9 = 28128954
pca.end.lrrc9 = 28169980

#### Pink lrrc9 #################################
pink_lrrc9_fst <- read.delim2("./results/fst/allele/pink-chum_NC_068455.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                              col.names = c("region", "chrName", "midPos", "Nsites", "Fst"),
                              row.names = NULL,sep = "\t") %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst) %>% ifelse(. < 0, 0, .), # remove negative Fst values
         chr = 35) %>%
  select(chr, midPos, Fst) %>%
  filter(midPos > xstart.lrrc9, midPos < xend.lrrc9) # cut Fst to start and end

#### Sockeye lrrc9 FST #####################
sock_lrrc9_fst <- read.delim("./results/fst/allele/sock-all_NC_068455.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                             col.names = c("region", "chrName", "midPos", "Nsites", "Fst"),
                             row.names = NULL,sep = "\t") %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst) %>% ifelse(. < 0, 0, .), # remove negative Fst values
         chr = 35) %>%
  select(chr, midPos, Fst) %>%
  filter(midPos > xstart.lrrc9, midPos < xend.lrrc9) # cut Fst to start and end

#### Chum lrrc9 ######################
chum_lrrc9_fst <- read.delim2("./results/fst/allele/chumrun_NC_068455.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                              col.names = c("region", "chrName", "midPos", "Nsites", "Fst"),
                              row.names = NULL,sep = "\t") %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst) %>% ifelse(. < 0, 0, .), # remove negative Fst values
         chr = 35) %>%
  select(chr, midPos, Fst) %>%
  filter(midPos > xstart.lrrc9, midPos < xend.lrrc9) # cut Fst to start and end

#### LRRC9 Genes ########################################################

# only the region of interest
gff_chr35_region <- gff_chr35 %>%
  mutate(start.pos = start.pos/1e6, fin.pos = fin.pos/1e6) %>% 
  filter(fin.pos > xstart.lrrc9, start.pos < xend.lrrc9)

# create new columns for genes and exons from ID
gff_exon_lrrc9 <- gff_chr35_region %>%
  filter(exon == "exon") %>%
  mutate(gene = str_match(ID, gene_pattern)[,2],  # gene abbr.      
         exonID = str_match(ID,exon_pattern)[,2]) # mRNA name and exon number

# this exon in six6a is too small that it doesn't even plot
# make slightly larger so it is visible in plot
gff_exon_lrrc9$fin.pos[which(gff_exon_lrrc9$gene == "six6a")[1]] <- 27.994500 # changed from 27.994277

exons_plot_chr35 <- gff_exon_lrrc9[,c(4,5,7:8)] # only retain columns of interest
# exons_plot_chr35$gene <- factor(exons_plot_chr35$gene, levels = unique(exons_plot_chr35$gene)) # factor based on gene name

# gff_chr35_region %>% 
#     filter(exon == "gene") %>%
#     mutate(gene = str_match(ID, gene_pattern)[,2],
#            geneName = str_match(ID, descr_pattern)[,2]) %>%
#     select(chrName, start.pos, fin.pos, gene, geneName) %>%
#     write_tsv("./data/R/chr35_chum_genes_all.tsv", col_names = T)

# The below file is manually edited for plotting purposes
# Colors are assigned to each gene
# genes_chr35_df <- read.csv("./data/R/chr35_chum_genes_exons_new.csv", header = T, 
#                            sep = "\t", row.names = NULL) %>%
#   mutate(beg.pos = as.numeric(beg.pos), end.pos = as.numeric(end.pos),
#          y.min = as.numeric(y.min), y.max = as.numeric(y.max))

genes_chr35_df <- read_tsv("./data/R/chr35_chum_genes_filtered.tsv", 
                           col_names = T, show_col_types = F) %>%
  mutate(beg.pos = as.numeric(beg.pos), end.pos = as.numeric(end.pos),
         y.min = as.numeric(y.min), y.max = as.numeric(y.max))

genes_chr35_df$geneAbbr <- factor(genes_chr35_df$geneAbbr, levels = genes_chr35_df$geneAbbr)
exons_plot_chr35 <- inner_join(gff_exon_lrrc9, genes_chr35_df %>% select(-chrName), by = "gene")
exons_plot_chr35$geneAbbr <- factor(exons_plot_chr35$geneAbbr, levels = unique(exons_plot_chr35$geneAbbr)) # factor based on gene name

# set factors for plotting columns
genes_chr35_df$gene <- factor(genes_chr35_df$gene, levels = genes_chr35_df$gene)

###### Make genes with FST < 0.5 gray w/o legend ##############
highfst_lrrc9 <- rbind(pink_lrrc9_fst, sock_lrrc9_fst, chum_lrrc9_fst) %>%
  filter(Fst > 0.5) %>%
  distinct(midPos) %>%
  pull(midPos)

lowfst_lrrc9_genes <- gff_chr35_region %>%
  filter(exon == 'gene') %>%
  mutate(gene = str_match(ID, gene_pattern)[,2]) %>%  # gene abbr.      
  rowwise() %>%
  mutate(highfst = any(highfst_lrrc9 >= start.pos & highfst_lrrc9 <= fin.pos)) %>%
  ungroup() %>%
  filter(highfst == F)

# don't include LOC### here
highfst_lrrc9_genes <- genes_chr35_df %>%
  filter(!(gene %in% lowfst_lrrc9_genes$gene) | grepl("^six",gene),
         !grepl("^LOC",gene))

# Include LOC### here
lowfst_lrrc9_genes <- genes_chr35_df %>%
  filter((gene %in% lowfst_lrrc9_genes$gene) | grepl("^LOC",gene),
         !grepl("^six",gene))

# only keep exons from genes that have color codes
highfst_lrrc9_exons <- filter(exons_plot_chr35, gene %in% highfst_lrrc9_genes$gene)
lowfst_lrrc9_exons <- filter(exons_plot_chr35, gene %in% lowfst_lrrc9_genes$gene)

##### Plotting ###########################

mypalette35 <- genes_chr35_df$color
names(mypalette35) <- levels(genes_chr35_df$gene)

# Set the general themes
theme_set(
  theme(legend.position = "none", 
    axis.text.x = element_text(size = smallFS, angle = 0, color = "black", vjust = 0.5),
    axis.title.x = element_text(size = mediumFS, angle = 0, color = "black"),
    axis.text.y = element_text(size = smallFS, angle = 0, color = "black", vjust = 0.5),
    axis.title.y = element_text(size = mediumFS, angle = 90,
                                margin = margin(t = 0, r = 4, b = 0, l = 0)),
    strip.text.y = element_text(angle = 0), 
    panel.grid.major = element_line(color = "gray95"),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(), 
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = "NA"),
    panel.background = element_rect(fill = "white")
  )
)


####### plot genes/exons
gene_lrrc9_plot <- ggplot() +
  # gray blocks without legend
  geom_rect(data = lowfst_lrrc9_genes, aes(xmin = beg.pos, xmax = end.pos, ymin = y.min, ymax = y.max),
            fill = "gray80", alpha = 0.7) +
  geom_rect(data = lowfst_lrrc9_exons, aes(xmin = start.pos, xmax = fin.pos, ymin = y.min-0.02, ymax = y.max+0.02),
            fill = "gray80", alpha = 0.7) +
  geom_rect(data = highfst_lrrc9_genes, aes(xmin = beg.pos, xmax = end.pos, ymin = y.min, ymax = y.max,
                                            fill = gene), alpha = 1) +
  geom_rect(data = highfst_lrrc9_exons, aes(xmin = start.pos, xmax = fin.pos, ymin = y.min-0.02, ymax = y.max+0.02,
                                            fill = gene)) +
  scale_y_continuous(limits = c(0, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.lrrc9, xend.lrrc9)) +
  scale_fill_manual(values = mypalette35)+   
  theme_minimal() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    strip.text.y = element_blank(), panel.spacing = unit(0.1,"lines"),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.text = element_text(size = smallFS, face = "italic"), 
    legend.title = element_text(size = mediumFS),
    # legend.key.size = unit(0.15, "in"),
    # legend.spacing.y = unit(0.0, "in"),
    # legend.justification = c(0, 1),
    axis.line = element_blank(), plot.margin = margin(0.0,0,0,0, "cm")) +
  guides(fill = guide_legend(title = "Genes",
                             keywidth  = unit(0.2, "in"),
                             keyheight = unit(0.15, "in")))
# gene_lrrc9_plot

# Pink lrrc9 Fst
pink_lrrc9_plot <- ggplot() +
  geom_vline(xintercept = pca.start.lrrc9/1e6, linetype = "dashed", linewidth = 0.8, color ="gray40") +
  geom_vline(xintercept = pca.end.lrrc9/1e6, linetype = "dashed", linewidth = 0.8, color ="gray40") +
  geom_point(data = pink_lrrc9_fst, aes(x = midPos, y = Fst), 
             size = pointSize, alpha = 0.6, color = "gray10") +
  labs(y = "Pink", x = "Chr 35 Position (Mb)") +
  scale_y_continuous(limits = c(-0.02, 1.02),
                     breaks = seq(0, 1, by = 0.5),
                     expand = expansion(mult = c(0.001, 0.01))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.lrrc9, xend.lrrc9),
                     breaks = seq(0, 100, by = 0.1)) +  
  theme(plot.margin = margin(0.1,0.15,0.1,0.05, "cm"))

# Sockeye lrrc9 Fst
sock_lrrc9_plot <- ggplot() +
  geom_vline(xintercept = pca.start.lrrc9/1e6, linetype = "dashed", linewidth = 0.8, color = "gray40") +
  geom_vline(xintercept = pca.end.lrrc9/1e6, linetype = "dashed", linewidth = 0.8, color = "gray40") +
  geom_point(data = sock_lrrc9_fst, aes(x = midPos, y = Fst), 
             size = pointSize, alpha = 0.6, color = "gray10") + 
  labs(y = "Sockeye") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(-0.02, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.lrrc9, xend.lrrc9),
                     breaks = seq(0, 100, by = 0.1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        plot.margin = margin(0.1,0.15,0.1,0.05, "cm"))

# Chum lrrc9 Fst
chum_lrrc9_plot <- ggplot() +
  geom_vline(xintercept = pca.start.lrrc9/1e6, linetype = "dashed", linewidth = 0.8, color = "gray40") +
  geom_vline(xintercept = pca.end.lrrc9/1e6, linetype = "dashed", linewidth = 0.8, color = "gray40") +
  geom_point(data = chum_lrrc9_fst, aes(x = midPos, y = Fst), 
             size = pointSize, alpha = 0.6, color = "gray10") +
  labs(y = "Chum") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(-0.01, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.lrrc9, xend.lrrc9),
                     breaks = seq(0, 100, by = 0.1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        plot.margin = margin(0.1,0.15,0.1,0.05, "cm"))

######### Combine Plots
multiplot_lrrc9_temp <- gene_lrrc9_plot / sock_lrrc9_plot / chum_lrrc9_plot / pink_lrrc9_plot + 
  plot_layout(heights = c(0.35, 1, 1, 1), guides = "collect") &
  theme(legend.position = 'right',
        legend.justification = 'top', legend.justification.right = c(0,1))

### 2B) ESRB FST ##############################################################

# which region of chr29 to plot (first and last position)
# panel spanning larger region (Fig 5)
xstart.esrb = 24.7
xend.esrb = 26.4
# where were boundaries for allele-based PCA
pca.start.esrb = 25414060
pca.end.esrb = 25501622

#### Chum esrb FST ##############################

chum_esrb_fst <- read.delim2("./results/fst/allele/chumrun_NC_068449.1_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                             col.names = c("region", "chrName", "midPos", "Nsites", "Fst"),
                             row.names = NULL,sep = "\t") %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst) %>% ifelse(. < 0, 0, .), # remove negative Fst values
         chr = 29) %>%
  select(chr, midPos, Fst) %>%
  filter(midPos > xstart.esrb, midPos < xend.esrb)

#### Coho esrb FST ###########################

coho_esrb_fst <- read.delim2("./results/fst/allele/coho-chum_NC_068449.1_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                             col.names = c("region", "chrName", "midPos", "Nsites", "Fst"),
                             row.names = NULL,sep = "\t") %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst) %>% ifelse(. < 0, 0, .), # remove negative Fst values
         chr = 29) %>%
  select(chr, midPos, Fst) %>%
  filter(midPos > xstart.esrb, midPos < xend.esrb) 

#### ESRB Genes ##############################################################

# only the region of interest
gff_chr29_region <- gff_chr29 %>%
  mutate(start.pos = start.pos/1e6, fin.pos = fin.pos/1e6) %>% 
  filter(fin.pos > xstart.esrb, start.pos < xend.esrb)

# create new columns for genes and exons from ID
gff_exon_esrb <- gff_chr29_region %>%
  filter(exon == "exon") %>%
  mutate(gene = str_match(ID, gene_pattern)[,2],  # gene abbr.      
         exonID = str_match(ID,exon_pattern)[,2]) # mRNA name and exon number

exons_esrb <- gff_exon_esrb[,c(4,5,7:8)] # only retain columns of interest

# only run this once
# genes_chr29_df <- gff_chr29_region %>%
#   filter(exon == "gene") %>%
#   mutate(gene = str_match(ID, gene_pattern)[,2],
#          geneName = str_match(ID, descr_pattern)[,2]) %>%
#   select(chrName, start.pos, fin.pos, gene, geneName) %>%
#   write_tsv("./data/R/chr29_chum_genes_all.tsv", col_names = T)

# this file was manually edited for plotting purposes based on above description
# Older Version: This filtered file removed all uncharacterized loci except those of great interest (esrb, six genes)
# genes_chr29_df <- read.csv("./data/R/chr29_chum_genes_exons_filtered.csv",
#                      header = T, row.names = NULL) %>%
#   mutate(beg.pos = as.numeric(beg.pos), end.pos = as.numeric(end.pos),
#          y.min = as.numeric(y.min), y.max = as.numeric(y.max))

# Newer Version with uncharacterized loci that will be grayed out
genes_chr29_df <- read_tsv("./data/R/chr29_chum_genes_filtered.tsv",
                           col_names = T, show_col_types = F)


genes_chr29_df$geneAbbr <- factor(genes_chr29_df$geneAbbr, levels = genes_chr29_df$geneAbbr) # set factors for plotting columns 
exons_plot_chr29 <- inner_join(exons_esrb, genes_chr29_df, by = "gene")
exons_plot_chr29$geneAbbr <- factor(exons_plot_chr29$geneAbbr, levels = unique(exons_plot_chr29$geneAbbr)) # factor based on gene name

unique(exons_plot_chr29$gene)

######### edit exons to plot - some are so small they don't register
exons_plot_chr29 <- exons_plot_chr29 %>% 
  mutate(beg.pos = ifelse(beg.pos < xstart.esrb, xstart.esrb, beg.pos),
         end.pos = ifelse(end.pos > xend.esrb, xend.esrb, end.pos),
         exon_length = (fin.pos - start.pos)*1e6,
         fin.pos = ifelse(exon_length > 1200, fin.pos, fin.pos + 0.5*(1200 - exon_length)/1e6),
         beg.pos = ifelse(exon_length > 1200, beg.pos, beg.pos - 0.5*(1200 - exon_length)/1e6),
         # run again
         fin.pos = ifelse(exon_length > 1200, fin.pos, fin.pos + 0.5*(1200 - exon_length)/1e6),
         beg.pos = ifelse(exon_length > 1200, beg.pos, beg.pos - 0.5*(1200 - exon_length)/1e6))

# this exon in esrb is too small even with the above edit
# make slightly larger so it is visible in plot
exons_plot_chr29$fin.pos[which(exons_plot_chr29$geneAbbr == "esrb")[length(which(exons_plot_chr29$geneAbbr == "esrb"))]] <- 25.456 # changed from 27.994277

##### Make genes with FST < 0.5 gray w/o legend ####################
highfst_esrb <- rbind(coho_esrb_fst, chum_esrb_fst) %>%
  filter(Fst > 0.5) %>%
  distinct(midPos)

lowfst_esrb_genes <- gff_chr29_region %>%
  filter(exon == 'gene') %>%
  mutate(gene = str_match(ID, gene_pattern)[,2]) %>%  # gene abbr.      
  rowwise() %>%
  mutate(highfst = any(highfst_esrb$midPos >= start.pos & highfst_esrb$midPos <= fin.pos)) %>%
  ungroup() %>%
  filter(highfst == F)

# Older Version
# highfst_esrb_genes <- genes_chr29_df %>% filter(!(gene %in% lowfst_esrb_genes$gene))
# lowfst_esrb_genes <- genes_chr29_df %>% filter((gene %in% lowfst_esrb_genes$gene))

# Newer Version
# don't include LOC### here
highfst_esrb_genes <- genes_chr29_df %>%
  filter(!(gene %in% lowfst_esrb_genes$gene),
         !grepl("^LOC",geneAbbr),
         !grepl("^si",geneAbbr))

# Include LOC### here
lowfst_esrb_genes <- genes_chr29_df %>%
  filter((gene %in% lowfst_esrb_genes$gene) | grepl("^LOC",geneAbbr) | grepl("^si",geneAbbr))


# only keep exons from genes that have color codes
highfst_esrb_exons <- filter(exons_plot_chr29, gene %in% highfst_esrb_genes$gene)
lowfst_esrb_exons <- filter(exons_plot_chr29, gene %in% lowfst_esrb_genes$gene)

##### Plotting #################

mypalette29 <- genes_chr29_df$color # color based on color column
names(mypalette29) <- levels(genes_chr29_df$geneAbbr)

# ESRB Gene plot
gene_esrb_plot <- ggplot() +
  geom_rect(data = lowfst_esrb_genes, aes(xmin = beg.pos, xmax = end.pos, ymin = y.min, ymax = y.max),
            fill = "gray80", alpha = 0.7) +
  geom_rect(data = lowfst_esrb_exons, aes(xmin = start.pos, xmax = fin.pos, ymin = y.min-0.02, ymax = y.max+0.02),
            fill = "gray80", alpha = 0.7) +
  geom_rect(data = highfst_esrb_genes, 
            aes(xmin = beg.pos, xmax = end.pos, ymin = y.min, ymax = y.max, fill = geneAbbr)) +
  geom_rect(data = highfst_esrb_exons, 
            aes(xmin = start.pos, xmax = fin.pos, ymin = y.min-0.02, ymax = y.max+0.02, fill = geneAbbr)) +
  scale_y_continuous(limits = c(0, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.esrb, xend.esrb)) +
  scale_fill_manual(values = mypalette29) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    strip.text.y = element_blank(), axis.line = element_blank(),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.spacing = unit(0.1,"lines"),
    legend.text = element_text(size = smallFS, face = "italic"), 
    legend.title = element_text(size = mediumFS),
    legend.position = "right", plot.margin = margin(0.05,0,0,0, "cm")) +
  guides(fill = guide_legend(title = "Genes",
         keywidth  = unit(0.2, "in"),
         keyheight = unit(0.15, "in"))) 

# Coho esrb plot
coho_esrb_plot <- ggplot() + 
  geom_vline(xintercept = pca.start.esrb/1e6, linetype = "dashed", linewidth = 0.8, color ="gray40") +
  geom_vline(xintercept = pca.end.esrb/1e6, linetype = "dashed", linewidth = 0.8, color ="gray40") +
  geom_point(data = coho_esrb_fst, aes(x = midPos, y = Fst), 
             size = pointSize, alpha = 0.6, color = "gray10") +
  labs(y = "Coho") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(-0.01, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.esrb, xend.esrb),
                     breaks = seq(0, 100, by = 0.2)) +
  theme(plot.margin = margin(0.1,0.15,0.05,0.05, "cm"),
        axis.title.x = element_blank(), axis.text.x = element_blank())

# Chum esrb plot
chum_esrb_plot <- ggplot() + 
  geom_vline(xintercept = pca.start.esrb/1e6, linetype = "dashed", linewidth = 0.8, color ="gray40") +
  geom_vline(xintercept = pca.end.esrb/1e6, linetype = "dashed", linewidth = 0.8, color ="gray40") +
  geom_point(data = chum_esrb_fst, aes(x = midPos, y = Fst), 
             size = pointSize, alpha = 0.6, color = "gray10") +
  labs(y = "Chum", x="Chr 29 Position (Mb)") +
  scale_y_continuous(limits = c(-0.02, 1.02),
                     breaks = seq(0, 1, by = 0.5),
                     expand = expansion(mult = c(0.001, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.esrb, xend.esrb),
                     breaks = seq(0, 100, by = 0.2)) +          
  theme(plot.margin = margin(0.15,0.05,0.1,0.05, "cm"))


# plot three figures on top of one another
multiplot_esrb_temp <- gene_esrb_plot / coho_esrb_plot / chum_esrb_plot + 
  plot_layout(heights = c(0.35, 1, 1),
              guides = "collect") & 
  theme(legend.position = "right", legend.justification = 'top')

### C) ESR1 FST ##############################################################

# which region of chr29 to plot (first and last position)
# panel spanning larger region (Fig 5)
xstart.er1 = 23
xend.er1 = 24.5

# where were boundaries for allele-based PCA
pca.start.er1 = 23264160
pca.end.er1 = 23645630

#### Sockeye er1 FST ##############################

sock_chr8_Fst <- read.delim2("./results/fst/allele/euclide_NC_068428.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                             col.names = c("region", "chrName", "midPos", "Nsites", "Fst"),
                             row.names = NULL,sep = "\t") %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst) %>% ifelse(. < 0, 0, .),
         chr = 8) %>%
  dplyr::select(chr, midPos, Fst)

# filter both pops to desired start and end point
sock_er1_fst <- sock_chr8_Fst %>%
  filter(midPos > xstart.er1, midPos < xend.er1)

#### Coho er1 FST ###########################

coho_chr8_Fst <- read.delim2("./results/fst/allele/coho-chum_NC_068428.1_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                             col.names = c("region", "chrName", "midPos", "Nsites", "Fst"),
                             row.names = NULL,sep = "\t") %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst) %>% ifelse(. < 0, 0, .),
         chr = 8) %>%
  dplyr::select(chr, midPos, Fst)

# filter both pops to desired start and end point
coho_er1_fst <- coho_chr8_Fst %>%
  filter(midPos > xstart.er1, midPos < xend.er1) 

# create new start/end pos based on Fst > 0.5 cutoff
xstart.er1 <- (filter(rbind(coho_er1_fst,sock_er1_fst), Fst > 0.5) %>% slice_min(midPos))[,"midPos"] - 0.05 # add buffer
xend.er1 <- (filter(rbind(coho_er1_fst,sock_er1_fst), Fst > 0.5) %>% slice_max(midPos))[,"midPos"] + 0.05 # add buffer

xstart.er1 = round(xstart.er1,digits=1)
xend.er1 = round(xend.er1,digits=1)

#### ER1 Genes ##############################################################

# only the region of interest
gff_chr8_region <- gff_chr8 %>%
  mutate(start.pos = start.pos/1e6, fin.pos = fin.pos/1e6) %>% 
  filter(fin.pos > xstart.er1, start.pos < xend.er1)

# create new columns for genes and exons from ID
exons_er1 <- gff_chr8_region %>%
  filter(exon == "exon") %>%
  mutate(gene = str_match(ID, gene_pattern)[,2],  # gene abbr.      
         exonID = str_match(ID,exon_pattern)[,2]) %>% # mRNA name and exon number
  select(c(4,5,7:8)) # only retain columns of interest

# only run this once
# genes_chr8_df <- gff_chr8_region %>% 
#   filter(exon == "gene") %>%
#   mutate(gene = str_match(ID, gene_pattern)[,2],
#          geneName = str_match(ID, descr_pattern)[,2]) %>%
#   select(chrName, start.pos, fin.pos, gene, geneName) %>%
#   write_tsv("./data/R/chr8_chum_genes_exons.tsv", col_names = T)

# this file was manually edited for plotting purposes based on above description
# This filtered file removed all uncharacterized loci except those of great interest (er1, six genes)
genes_chr8_df <- read_tsv("./data/R/chr8_chum_genes_exons_filtered.tsv",
                          col_names =T,show_col_types = F) %>%
  mutate(beg.pos = as.numeric(beg.pos), end.pos = as.numeric(end.pos),
         y.min = as.numeric(y.min), y.max = as.numeric(y.max))

genes_chr8_df$geneAbbr <- factor(genes_chr8_df$geneAbbr, levels = unique(genes_chr8_df$geneAbbr)) # set factors for plotting columns 
exons_plot_chr8 <- inner_join(exons_er1, genes_chr8_df, by = "gene")
exons_plot_chr8$geneAbbr <- factor(exons_plot_chr8$geneAbbr, levels = unique(exons_plot_chr8$geneAbbr)) # factor based on gene name

unique(exons_plot_chr8$gene)

######### edit exons to plot - some are so small they don't register
exons_plot_chr8 <- exons_plot_chr8 %>% 
  mutate(beg.pos = ifelse(beg.pos < xstart.er1, xstart.er1, beg.pos),
         end.pos = ifelse(end.pos > xend.er1, xend.er1, end.pos),
         exon_length = (fin.pos - start.pos)*1e6,
         fin.pos = ifelse(exon_length > 1200, fin.pos, fin.pos + 0.5*(1200 - exon_length)/1e6),
         beg.pos = ifelse(exon_length > 1200, beg.pos, beg.pos - 0.5*(1200 - exon_length)/1e6),
         # run again
         fin.pos = ifelse(exon_length > 1200, fin.pos, fin.pos + 0.5*(1200 - exon_length)/1e6),
         beg.pos = ifelse(exon_length > 1200, beg.pos, beg.pos - 0.5*(1200 - exon_length)/1e6))

##### Make genes with FST < 0.5 gray w/o legend ####################
highfst_er1 <- rbind(coho_er1_fst, sock_er1_fst) %>%
  filter(Fst > 0.5) %>%
  distinct(midPos)

lowfst_er1_genes <- gff_chr8_region %>%
  filter(exon == 'gene') %>%
  mutate(gene = str_match(ID, gene_pattern)[,2]) %>%  # gene abbr.      
  rowwise() %>%
  mutate(highfst = any(highfst_er1$midPos >= start.pos & highfst_er1$midPos <= fin.pos)) %>%
  ungroup() %>%
  filter(highfst == F)

# don't include LOC### or si:## genes here
highfst_er1_genes <- genes_chr8_df %>%
  filter(!(gene %in% lowfst_er1_genes$gene),
         !grepl("^LOC",geneAbbr),
         !grepl("^si",geneAbbr))

# Include LOC### and si:## genes here
lowfst_er1_genes <- genes_chr8_df %>%
  filter((gene %in% lowfst_er1_genes$gene) | grepl("^LOC",geneAbbr) | grepl("^si",geneAbbr))

# only keep exons from genes that have color codes
highfst_er1_exons <- filter(exons_plot_chr8, gene %in% highfst_er1_genes$gene)
lowfst_er1_exons <- filter(exons_plot_chr8, gene %in% lowfst_er1_genes$gene)

#### Plotting #################

mypalette08 <- genes_chr8_df$color # color based on color column
names(mypalette08) <- levels(genes_chr8_df$geneAbbr)

# Genes/exons
er1_gene_plot <- ggplot() +
  geom_rect(data = lowfst_er1_genes, aes(xmin = beg.pos, xmax = end.pos, ymin = y.min, ymax = y.max),
            fill = "gray80", alpha = 0.7) +
  geom_rect(data = lowfst_er1_exons, aes(xmin = start.pos, xmax = fin.pos, ymin = y.min-0.02, ymax = y.max+0.02),
            fill = "gray80", alpha = 0.7) +
  geom_rect(data = highfst_er1_genes, aes(xmin = beg.pos, xmax = end.pos, ymin = y.min, ymax = y.max,
                                          fill = geneAbbr)) +
  geom_rect(data = highfst_er1_exons, aes(xmin = start.pos, xmax = fin.pos,ymin = y.min-0.02, ymax = y.max+0.02,
                                          fill = geneAbbr)) +
  scale_y_continuous(limits = c(0, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.er1, xend.er1)) +
  scale_fill_manual(values = mypalette08) +
  theme_minimal() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    strip.text.y = element_blank(), axis.line = element_blank(),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.spacing = unit(0.1,"lines"),
    legend.text = element_text(size = smallFS, face = "italic"), 
    legend.title = element_text(size = mediumFS), 
    legend.position = "right", plot.margin = margin(0.05,0,0,0, "cm")) +
  guides(fill = guide_legend(title = "Genes",
                             keywidth  = unit(0.2, "in"),
                             keyheight = unit(0.15, "in")))   

# Sockeye ESR1 FST
sock_er1_plot <- ggplot() + 
  geom_vline(xintercept = pca.start.er1/1e6, linetype = "dashed", linewidth = 0.8, color = "gray40") +
  geom_vline(xintercept = pca.end.er1/1e6, linetype = "dashed", linewidth = 0.8, color = "gray40") +
  geom_point(data = sock_er1_fst, aes(x = midPos, y = Fst), 
             size = pointSize, alpha = 0.6, color = "gray10") +
  labs(y = "Sockeye") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(-0.01, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.er1, xend.er1),
                     breaks = seq(0, 100, by = 0.2)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        plot.margin = margin(0.1,0.15,0.05,0.05, "cm"))

# Coho ESR1 FST 
coho_er1_plot <- ggplot() +
  geom_vline(xintercept = pca.start.er1/1e6, linetype = "dashed", linewidth = 0.8, color = "gray40") +
  geom_vline(xintercept = pca.end.er1/1e6, linetype = "dashed", linewidth = 0.8, color = "gray40") +
  geom_point(data = coho_er1_fst, aes(x = midPos, y = Fst), 
             size = pointSize, alpha = 0.6, color = "gray10") +
  labs(x="Chr 8 Position (Mb)", y = "Coho") +
  scale_y_continuous(limits = c(-0.02, 1.02),
                     breaks = seq(0, 1, by = 0.5),
                     expand = expansion(mult = c(0.001, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.er1, xend.er1),
                     breaks = seq(0, 100, by = 0.2)) +          
  theme(plot.margin = margin(0.15,0.05,0.1,0.05, "cm"))

# plot three figures on top of one another
er1_multiplot_temp <- er1_gene_plot / sock_er1_plot / coho_er1_plot + 
  plot_layout(heights = c(0.35, 1, 1),
              guides = "collect") & 
  theme(legend.position = "right", legend.justification = 'top')

### A/B/C  #####################

# add FST as separate label (it will be it's own plot)
y_lab <- ggplot() +
  annotate(geom = "text", x = 1, y = 1, angle = 90,
           label = expression(italic(F[ST])),
           size = largeFS / .pt) + # need to convert annotation to pt 
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0,0.05,0,0,unit="in")) 

# combine FST label to other plots 
multiplot_lrrc9 <- (y_lab - multiplot_lrrc9_temp) + # patchwork uses hyphen to allow for lefthand additions
  plot_layout(widths = c(1, 40))

multiplot_esrb <- (y_lab - multiplot_esrb_temp) + # patchwork hyphen adds to left
  plot_layout(widths = c(1, 40))

multiplot_er1 <- (y_lab - er1_multiplot_temp) + # patchwork uses hyphen to allow for lefthand additions
  plot_layout(widths = c(1, 40))

# # for viewing
# multiplot_lrrc9
# multiplot_esrb
# multiplot_er1

esr1_title <- ggdraw() + 
  draw_label(expression(bold(A)~~~italic(esr1)),
             x = 0, hjust = 0, vjust=0.75,
             fontfamily = "sans",
             size = largeFS) + 
  theme(plot.margin = margin(0, 0, 0, 0))

esrb_title <- ggdraw() + 
  draw_label(expression(bold(B)~~~italic(esr2b)),
             x = 0, hjust = 0, vjust=0.75,
             fontfamily = "sans",
             size = largeFS) + 
  theme(plot.margin = margin(0, 0, 0, 0))

lrrc9_title <- ggdraw() + 
  draw_label(expression(bold(C)~~~italic(lrrc9)),
             x = 0, hjust = 0, vjust=0.75,
             fontfamily = "sans",
             size = largeFS) + 
  theme(plot.margin = margin(0, 0, 0, 0))

threefst_again <- (esr1_title /multiplot_er1 / esrb_title / multiplot_esrb / lrrc9_title / multiplot_lrrc9) +
  plot_layout(heights = c(0.2, 2, 0.2,2, 0.2, 2.9), ncol = 1)


png(paste0("./figures/fst/Figure5_Allele_FST_Genes_",format(Sys.Date(),"%Y%m%d"),".png"), 
     width = width2col, height = 1.3*width2col, res = 500, units = "in")
print(threefst_again)
dev.off()

# svg for science 
svglite(str_c("./figures/fst/Figure5_Allele_FST_Genes_",format(Sys.Date(),"%Y%m%d"),".svg"),
        width = width2col, height = 1.3*width2col)
print(threefst_again)
dev.off()
