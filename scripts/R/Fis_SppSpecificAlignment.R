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
  spp = c("Chum", "Chum", "Chum", "Pink", "Sockeye", 
          "Sockeye", "Coho", "Coho"),
  gene = c("esr1", "lrrc9", "esr2b", "lrrc9", "lrrc9", 
           "esr1", "esr1", "esr2b"),
  pos = c("8:23523137-23555585", "35:28124851-28174083", "29:25448834-25486645", "10:72708379-72737558", "12:41182696-41210041",
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

all_chrom_meta <- bind_rows(chrom_dfs, .id = "spp")

gene_by_spp <- genes_df %>%
  mutate(chr = as.character(chr)) %>%
  left_join(all_chrom_meta, join_by('spp','chr'))

# bring in fis data by chromosomes of interest as defined by genes above
fis_df <- map2_dfr(gene_by_spp$spp, gene_by_spp$chrName, function(species, chr) {
  
  file_paths <- list.files("./results/fis", 
                           pattern = str_c(".*", chr, ".*HWE_minMaf\\.hwe\\.gz$"),
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

# should be 2 of each for early/late
fis_df %>%
  distinct(spp,runtime,chrName) %>%
  count(spp,chrName)

### The plots are a lot, so maybe just do a table ##########################
boundaries <- read_tsv('./data/R/gene_peak_boundaries.txt', col_names=TRUE) %>%
  separate(Boundaries, into=c('start','end'),sep='-', convert = TRUE)

# chr 8 in chum (esr1) is included because it is in the mappability calcs.

fis_by_bound <- boundaries %>%
  mutate(Chromosome = as.character(Chromosome)) %>%
  left_join(fis_df %>%
              left_join(all_chrom_meta, join_by('spp','chrName')),
            by = join_by(Species == spp, Chromosome == chr, 
                         start <= Position, 
                         end >= Position)) %>%
  mutate(lowFis = Fis < -0.5) %>%
  group_by(Species, Chromosome, start, end) %>%
  summarize(Average = mean(Fis, na.rm = TRUE),
            Median = median(Fis, na.rm = TRUE),
            SD = sd(Fis, na.rm = TRUE),
            Min = min(Fis, na.rm = TRUE),
            Max = max(Fis, na.rm = TRUE),
            PropLowFis = mean(lowFis, na.rm = TRUE), # mean is prop for T/F
            .groups = "drop")

boundaries %>%
  mutate(Chromosome = as.character(Chromosome)) %>%
  left_join(fis_by_bound, join_by(Species, Chromosome, start, end)) %>%
  write_tsv("./results/TableSX_FIS_by_gene_region.tsv")

fis_df %>%
  filter(spp=='Coho',
         chr=='21',
         Position >= 24719119, Position <= 24743990)


##### Plotting ###########################

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
    
    ggsave(filename = str_c("./figures/fis/minMaf",spp,"_chr",chr,"_",gene,".jpeg"),
           plot = p, width = 7, height = 5, dpi = 200)
  }
)


### looking at the data
boundaries %>%
  mutate(Chromosome = as.character(Chromosome)) %>%
  left_join(fis_df %>%
              left_join(all_chrom_meta, join_by('spp','chrName')),
            by = join_by(Species == spp, Chromosome == chr, 
                         start <= Position, end >= Position)) %>% View
