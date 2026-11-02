# Species-specific Alignment FIS 

packages_needed <- c("ggplot2", "scales", "ggpubr", "ggrepel", "stringr", 
                     "data.table", "plyr","tools","gtools","reshape2", 
                     "patchwork", "cowplot", "tidyverse")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}
rm(i, packages_needed)

# manually create data frame with information on three genes with shared peaks and their gene boundaries
genes_df <- data.frame(
  spp = c("Chum", "Chum", "Pink", "Sockeye", "Sockeye", "Coho", "Coho"),
  gene = c("lrrc9", "esr2b", "lrrc9", "lrrc9", "esr1", "esr1", "esr2b"),
  pos = c("35:28124851-28174083", "29:25448834-25486645", "10:72708379-72737558", "12:41182696-41210041",
          "13:7619247-7647640", "21:24716632-24746477", "7:30261985-30294176")
) %>%
  separate(pos, into = c("chr", "start", "end"), 
           sep = "[:\\-]", convert = TRUE)

chrfiles <- c(Chum = "./data/R/chumV2_chrom_meta.txt", # no header
  Pink = "./data/R/pinkV1_chrom_meta.txt",
  Sockeye = "./data/R/sockeyeV1_chrom_meta.txt",
  Coho = "./data/R/cohoV2_chrom_meta.txt")

# Read and rename columns in one step
chrom_dfs <- lapply(chrfiles, function(file) {
  read.table(file, header = FALSE, col.names = c("chr", "chrName")) %>%
    mutate(chr = as.character(chr))
})

localscorefiles <- c(Chum = "./results/localscore/chumrun_minInd0.3_minDepthHalf_MM3_xi3_01sig.txt",
                     Pink = "./results/localscore/AukePink_minInd0.3_HWE_xi3_01sig.txt",
                     Sockeye = "./results/localscore/sock-sock_minInd0.3_HWE_xi3_01sig.txt",
                     Coho = "./results/localscore/tocoho_minInd0.3_minDepthHalf_HWE_xi3_01sig.txt")

ls_peaks <- lapply(localscorefiles, function(file) {
  read_tsv(file, col_names = T, show_col_types=F) %>%
  dplyr::rename(`chrName` = `chr`)
})

# add normal chr column by species
spp_spec_ls_peaks <- Map(function(ls, chrom) {
  left_join(ls, chrom, by = "chrName")
}, ls_peaks, chrom_dfs)

all_chrom_meta <- bind_rows(chrom_dfs, .id = "spp")

gene_by_spp <- genes_df %>%
  mutate(chr = as.character(chr)) %>%
  left_join(all_chrom_meta, join_by('spp','chr'))

# bring in fis data by chromosomes of interest as defined by genes above
fis_df <- map2_dfr(gene_by_spp$spp, gene_by_spp$chrName, function(species, chr) {
  
  file_paths <- list.files("./results/fis", 
                           pattern = str_c(".*", chr, ".*\\.hwe\\.gz$"),
                           full.names = TRUE)
  
  map_dfr(file_paths, function(fp) {
    rt <- basename(fp) %>% str_split_i("_", 4)
    
    read_tsv(fp, col_names = TRUE, show_col_types = FALSE) %>%
      mutate(spp = species,
             runtime = rt,
             midPos = Position/1e6,
             Fis = as.numeric(`F`),
             log10p  = -log10(as.numeric(`p-value`))) %>%
      dplyr::select(-c(LRT, `F`, hweFreq, Freq, `p-value`)) %>%
      dplyr::rename(chrName = Chromo)
  })
})

fis_df %>%
  distinct(spp,runtime,chrName)

##### Plotting ###########################

plot_list <- pmap(gene_by_spp, 
                  function(spp, chrName, chr, start, end, gene) {
    fis_df %>%
      filter(spp == !!spp, chrName == !!chrName) %>%
      ggplot() +
      annotate("rect", xmin=start/1e6,xmax=end/1e6,ymin=-Inf,ymax=Inf,
                alpha=0.5,fill="khaki1")+
      geom_point(aes(x = midPos, y = Fis), alpha = 0.8) +
      coord_cartesian(xlim = c(start/1e6-0.5, end/1e6+0.5)) + # in Mb
      facet_grid(runtime~.) +
      labs(title = str_c(spp, " (", gene, ")"),
           x = paste("Chromosome", chr, "Position (Mb)"),
           y = expression(italic(F)[IS])) +
      theme_bw()
  }
)

pwalk(
  gene_by_spp,
  function(spp, chrName, chr, start, end, gene, ...) {
    p <- fis_df %>%
      filter(spp == !!spp, 
             chrName == !!chrName) %>%
      ggplot() +
      annotate("rect", xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf,
               alpha = 0.5, fill = "khaki1") +
      geom_point(aes(x = midPos, y = Fis), alpha = 0.8) +
      coord_cartesian(xlim = c(start / 1e6 - 0.5, end / 1e6 + 0.5)) +
      facet_grid(runtime ~ .) +
      labs(title = str_c(spp, " (", gene, ")"),
           x = paste("Chromosome", chr, "Position (Mb)"),
           y = expression(italic(F)[IS])) +
      theme_bw()
    
    ggsave(filename = str_c("./figures/fis/",spp,"_chr",chr,"_",gene,".jpeg"),
           plot = p, width = 7, height = 5, dpi = 500)
  }
)





####### plot genes/exons
gene_lrrc9_plot <- ggplot() +
  # gray blocks without legend
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.lrrc9, xend.lrrc9)) +
  scale_fill_manual(values = mypalette) +
  guides(fill = guide_legend(title = "Genes")) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    strip.text.y = element_blank(), panel.spacing = unit(0.1,"lines"),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.text = element_text(size = 24), legend.title = element_text(size = 26),
    axis.line = element_blank(), plot.margin = unit(c(0.05,0,0,0), "cm"))
# gene_lrrc9_plot

# Pink lrrc9 Fst
pink_lrrc9_plot <- ggplot() +
  geom_vline(xintercept = pca.start.lrrc9/1e6, linetype = "dashed", linewidth=1.2) +
  geom_vline(xintercept = pca.end.lrrc9/1e6, linetype = "dashed", linewidth=1.2) +
  geom_point(data = pink_lrrc9_fst, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") +
  labs(y = "Pink", x = "Chr 35 Position (Mb)") +
  scale_y_continuous(limits = c(-0.02, 1.02),
                     breaks = seq(0, 1, by = 0.5),
                     expand = expansion(mult = c(0.001, 0.01))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.lrrc9, xend.lrrc9),
                     breaks = seq(0, 100, by = 0.1)) +  
  theme(#axis.text.x = element_text(angle = 0, size = 18, color = "black", vjust = 0.5),
        #axis.title.x = element_text(size = 20), strip.text.y = element_text(angle = 0),
        plot.margin = unit(c(0.1,0.15,0.1,0.05), "cm"))

# Sockeye lrrc9 Fst
sock_lrrc9_plot <- ggplot() +
  geom_vline(xintercept = pca.start.lrrc9/1e6, linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = pca.end.lrrc9/1e6, linetype = "dashed", linewidth = 1.2) +
  geom_point(data = sock_lrrc9_fst, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") + 
  labs(y = "Sockeye") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(-0.02, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.lrrc9, xend.lrrc9),
                     breaks = seq(0, 100, by = 0.1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        plot.margin = unit(c(0.1,0.15,0.1,0.05), "cm"))

# Chum lrrc9 Fst
chum_lrrc9_plot <- ggplot() +
  geom_vline(xintercept = pca.start.lrrc9/1e6, linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = pca.end.lrrc9/1e6, linetype = "dashed", linewidth = 1.2) +
  geom_point(data = chum_lrrc9_fst, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") +
  labs(y = "Chum") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(-0.01, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.lrrc9, xend.lrrc9),
                     breaks = seq(0, 100, by = 0.1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        plot.margin = unit(c(0.1,0.15,0.1,0.05), "cm"))

######### Combine Plots
multiplot_lrrc9_temp <- gene_lrrc9_plot / sock_lrrc9_plot / chum_lrrc9_plot / pink_lrrc9_plot + 
  plot_layout(heights = c(0.35, 1, 1, 1), guides = "collect") &
  theme(legend.position = 'right',
        legend.justification = 'top', legend.justification.right = c(0,0.8))

# add FST as separate label (it will be it's own plot)
y_lab <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = expression(italic(F[ST])), angle = 90, size = 10) + 
  coord_cartesian(clip = "off") +
  theme_void() 

# combine FST label to other plots 
multiplot_lrrc9 <- (y_lab - multiplot_lrrc9_temp) + # patchwork uses hyphen to allow for lefthand additions
  plot_layout(widths = c(1, 12))
multiplot_lrrc9


### 2B) ESRB FST ##############################################################

# which region of chr29 to plot (first and last position)
# panel spanning larger region (Fig 5)
xstart.esrb = 24.7
xend.esrb = 26.4
# where were boundaries for allele-based PCA
pca.start.esrb = 25414060
pca.end.esrb = 25501622

#### Chum esrb FST ##############################

chum_esrb_fst <- read.delim2("./results/fis/chumrun_NC_068449.1_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                        col.names = c("region", "chrName", "midPos", "Nsites", "Fst"),
                        row.names = NULL,sep = "\t") %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst) %>% ifelse(. < 0, 0, .), # remove negative Fst values
         chr = 29) %>%
  select(chr, midPos, Fst) %>%
  filter(midPos > xstart.esrb, midPos < xend.esrb)

#### Coho esrb FST ###########################

coho_esrb_fst <- read.delim2("./results/fiscoho-alignment/tocoho_NC_068449.1_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
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

mypalette <- genes_chr29_df$color # color based on color column
  names(mypalette) <- levels(genes_chr29_df$geneAbbr)

##### Plotting #################

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
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.esrb, xend.esrb)) +
  scale_fill_manual(values = mypalette) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    strip.text.y = element_blank(), axis.line = element_blank(),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.spacing = unit(0.1,"lines"),
    legend.text = element_text(size = 24), legend.title = element_text(size = 26),
    legend.position = "right", plot.margin = unit(c(0.05,0,0,0), "cm")) +
  guides(fill = guide_legend(title = "Genes")) 

  
# Coho esrb plot
coho_esrb_plot <- ggplot() + 
  geom_vline(xintercept = pca.start.esrb/1e6, linetype = "dashed", linewidth=1.2) +
  geom_vline(xintercept = pca.end.esrb/1e6, linetype = "dashed", linewidth=1.2) +
  geom_point(data = coho_esrb_fst, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") +
  labs(y = "Coho") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(-0.01, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.esrb, xend.esrb),
                     breaks = seq(0, 100, by = 0.2)) +
  theme(plot.margin = unit(c(0.1,0.15,0.05,0.05), "cm"),
        axis.title.x = element_blank(), axis.text.x = element_blank())

# Chum esrb plot
chum_esrb_plot <- ggplot() + 
  geom_vline(xintercept = pca.start.esrb/1e6, linetype = "dashed", linewidth=1.2) +
  geom_vline(xintercept = pca.end.esrb/1e6, linetype = "dashed", linewidth=1.2) +
  geom_point(data = chum_esrb_fst, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") +
  labs(y = "Chum", x="Chr 29 Position (Mb)") +
  scale_y_continuous(limits = c(-0.02, 1.02),
                     breaks = seq(0, 1, by = 0.5),
                     expand = expansion(mult = c(0.001, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.esrb, xend.esrb),
                     breaks = seq(0, 100, by = 0.2)) +          
  theme(plot.margin = unit(c(0.15,0.05,0.1,0.05), "cm"))


# plot three figures on top of one another
multiplot_esrb_temp <- gene_esrb_plot / coho_esrb_plot / chum_esrb_plot + 
  plot_layout(heights = c(0.32, 1, 1),
              guides = "collect") & 
  theme(legend.position = "right", legend.justification = 'top')

# combine FST label to other plots 
multiplot_esrb <- (y_lab - multiplot_esrb_temp) + # patchwork hyphen adds to left
  plot_layout(widths = c(1, 12))
multiplot_esrb


### C) ESR1 FST ##############################################################

# which region of chr29 to plot (first and last position)
# panel spanning larger region (Fig 5)
xstart.er1 = 23
xend.er1 = 24.5

# where were boundaries for allele-based PCA
pca.start.er1 = 23264160
pca.end.er1 = 23645630

#### Sockeye er1 FST ##############################

sock_chr8_Fst <- read.delim2("./results/fis/euclide_NC_068428.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
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

coho_chr8_Fst <- read.delim2("./results/fis/coho-alignment/tocoho_NC_068428.1_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
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

mypalette <- genes_chr8_df$color # color based on color column
  names(mypalette) <- levels(genes_chr8_df$geneAbbr)

#### Plotting #################

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
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.er1, xend.er1)) +
  scale_fill_manual(values = mypalette) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    strip.text.y = element_blank(), axis.line = element_blank(),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.spacing = unit(0.1,"lines"),
    legend.text = element_text(size = 24), legend.title = element_text(size = 26),
    legend.position = "right", plot.margin = unit(c(0.05,0,0,0), "cm")) +
  guides(fill = guide_legend(title = "Genes")) 
  
# Sockeye ESR1 FST
sock_er1_plot <- ggplot() + 
  geom_vline(xintercept = pca.start.er1/1e6, linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = pca.end.er1/1e6, linetype = "dashed", linewidth = 1.2) +
  geom_point(data = sock_er1_fst, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") +
  labs(y = "Sockeye") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(-0.01, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.er1, xend.er1),
                     breaks = seq(0, 100, by = 0.2)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.15,0.05,0.05), "cm"))

# Coho ESR1 FST 
coho_er1_plot <- ggplot() +
  geom_vline(xintercept = pca.start.er1/1e6, linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = pca.end.er1/1e6, linetype = "dashed", linewidth = 1.2) +
  geom_point(data = coho_er1_fst, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") +
  labs(x="Chr 8 Position (Mb)", y = "Coho") +
  scale_y_continuous(limits = c(-0.02, 1.02),
                     breaks = seq(0, 1, by = 0.5),
                     expand = expansion(mult = c(0.001, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.er1, xend.er1),
                     breaks = seq(0, 100, by = 0.2)) +          
  theme(plot.margin = unit(c(0.15,0.05,0.1,0.05), "cm"))

# plot three figures on top of one another
multiplot_temp <- er1_gene_plot / sock_er1_plot / coho_er1_plot + 
  plot_layout(heights = c(0.32, 1, 1),
              guides = "collect") & 
  theme(legend.position = "right", legend.justification = 'top')

# combine FST label to other plots 
multiplot_er1 <- (y_lab - multiplot_temp) + # patchwork uses hyphen to allow for lefthand additions
  plot_layout(widths = c(1, 12))
multiplot_er1


### A/B/C  #####################
threefst_cowplot <- plot_grid(NULL, multiplot_er1, NULL, multiplot_esrb, NULL, multiplot_lrrc9,
                            rel_heights = c(0.2,2,0.2,2,0.2,2.8), ncol = 1, nrow = 6, 
                            align = 'v', axis = "lr", # axis is what allowed actually plot alignment
                            labels = c('A) esr1','','B) esr2b','','C) lrrc9',''), label_fontfamily = "helvetica",
                            label_size = 30, label_colour = "black")
# threefst_cowplot
jpeg(paste0("../Salmon_runtiming/2024_fourspecies/figures/Figure4a-c_FST_Allele_Genes_",format(Sys.Date(),"%Y%m%d"),".jpg"), 
#jpeg(paste0("./figures/Figure4a-c_Allele_FST_Genes_",format(Sys.Date(),"%Y%m%d"),".jpg"), 
          width = 19, height = 22, res = 500, units = "in")
print(threefst_cowplot)
dev.off()

### D/E/F Instead ##################

threefst_cowplot <- plot_grid(NULL, multiplot_er1, NULL, multiplot_esrb, NULL, multiplot_lrrc9,
                              rel_heights = c(0.2,2,0.2,2,0.2,2.8), ncol = 1, nrow = 6, 
                              align = 'v', axis = "lr", # axis is what allowed actually plot alignment
                              labels = c('D) esr1','','E) esr2b','','F) lrrc9',''), label_fontfamily = "helvetica",
                              label_size = 30, label_colour = "black")
# threefst_cowplot
jpeg(paste0("../Salmon_runtiming/2024_fourspecies/figures/Figure3d-f_FST_Allele_Genes_",format(Sys.Date(),"%Y%m%d"),".jpg"), 
#jpeg(paste0("./figures/Figure4a-c_Allele_FST_Genes_",format(Sys.Date(),"%Y%m%d"),".jpg"), 
     width = 19, height = 22, res = 500, units = "in")
print(threefst_cowplot)
dev.off()

