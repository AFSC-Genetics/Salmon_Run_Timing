# plot WHOLE GENOME pink even, sockeye, chum, and coho

packages_needed <- c("ggplot2", "scales", "ggpubr", "ggrepel", "tidyverse",
                     "stringr", "data.table", "plyr","gtools","reshape2", 
                     "here", "magrittr", "patchwork", "cowplot")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}
rm(packages_needed, i)  

# Table for chromosome name and number
chrom_df <- read.table("./data/R/chrom_meta.txt", header = TRUE)

###### Read in each per-chrom Fst File #########################################

########## PINK
pink_list <- as.list(Sys.glob("./results/fst/pink/pink-chum*early-late*minInd0.3*fst.txt"))

# read in all data files that match wildcard
pink_df <- pink_list %>%
  set_names(nm = pink_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL, 
                  col.names = c("region", "chrName", "midPos", "Nsites", "Fst"))
  ) %>% 
  select(-c(region, Nsites)) %>%
  mutate(Species = "Pink Even")

########## SOCKEYE
sock_list <- as.list(Sys.glob("./results/fst/sockeye/*euclide*minInd0.3*fst.txt"))

# read in all data files that match wildcard
sock_df <- sock_list %>%
  set_names(nm = sock_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL, 
                  col.names = c("region", "chrName", "midPos", "Nsites", "Fst"))
  ) %>% 
  select(-c(region, Nsites)) %>%
  mutate(Species = "Sockeye")

########## CHUM
chum_list <- as.list(Sys.glob("./results/fst/chum/*chumrun*minInd0.3_minDepth*fst.txt"))

# read in all data files that match wildcard
chum_df <- chum_list %>%
  set_names(nm = chum_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL, 
                  col.names = c("region", "chrName", "midPos", "Nsites", "Fst"))
  ) %>% 
  select(-c(region, Nsites)) %>%
  mutate(Species = "Chum")

########## COHO
coho_list <- as.list(Sys.glob("./results/fst/coho/*coho-chum*minInd0.3_minDepth*fst.txt"))

# read in all data files that match wildcard
coho_df <- coho_list %>%
  set_names(nm = coho_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL,
                  col.names = c("region", "chrName", "midPos", "Nsites", "Fst"))
  ) %>% 
  select(-c(region, Nsites)) %>%
  mutate(Species = "Coho")

############## COMBINE FOUR SPECIES INTO ONE DATAFRAME #######################
# bind rows together
four_df <- bind_rows(pink_df, sock_df, chum_df, coho_df)

# Change negative Fst SNPs to 0s
four_df$Fst[four_df$Fst < 0] <- 0

four_df <- four_df %>%
  mutate(midPos = as.numeric(midPos/10^6),
         Fst = as.numeric(Fst))

#add chromosome numbers, remove chrName
four_df <- left_join(four_df, chrom_df, by = "chrName") %>%
  select(chr, midPos, Fst, Species)
head(four_df)

rm(pink_df, sock_df, chum_df, coho_df, # remove to save space
   pink_list, sock_list, chum_list, coho_list) 

# # find max of each chrom
cumulate <- four_df %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(max_bp = max(midPos), .groups = "drop") %>%
  mutate(bp_offset = lag(cumsum(max_bp + 5), default = 0))

cum_df <- four_df %>%
  left_join(cumulate %>%
              dplyr::select(chr, bp_offset), by = "chr") %>%
  mutate(cumPos = midPos + bp_offset) %>%
  dplyr::select(-bp_offset)

rm(four_df) # remove unneeded files

# SINGLE OUT ESR1, ESRB, AND LRRC9 REGIONS
esr1_df <- cum_df %>%
  filter(chr == 8,
         midPos >= 23.141735,
         midPos <= 23.906114)

esrb_df <- cum_df %>%
  filter(chr == 29,
         midPos >= 25,
         midPos <= 26.2)

lrrc9_df <- cum_df %>%
  filter(chr == 35,
         midPos >= 27.85,
         midPos <= 28.20)

############# PREPARE FOR PLOTTING #########################################

axis_set <- cum_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos)) %>%
  filter(chr %% 2 == 1)

# Bc plotting on a massive scale, the loss of FST==0 will not change the figure.
cum_df <- cum_df %>%
  filter(Fst > 0)

# test pruning fst values below the threshold where it is just overlapping points
# this was tested across different quantiles and these do not change the appearance of the plot
cum_df <- cum_df %>%
  dplyr::group_by(Species) %>%
  filter(Fst >= quantile(Fst, 0.95, na.rm = TRUE) |
           (Fst < quantile(Fst, 0.95, na.rm = TRUE) & runif(n()) < 0.2)
  ) %>%
  ungroup()

# color blind palette colors
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999")

# for science, consistent font sizing
width2col <- 7.25  # 2 col, in inches

largeFS <- 10 # Panel Labels
mediumFS <- 8  # Axis titles, legend titles
smallFS <- 6  # Axis text, annotations, captions

# Set the ggplot theme
theme_set(
  theme( 
    legend.position = "none",
    panel.grid.major.y = element_blank(), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = smallFS, angle = 0),
    axis.title.x = element_text(size = mediumFS, 
                                margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = mediumFS, angle = 90, 
                                margin = margin(t = 0, r = 5, b = 0, l = 0)),
    panel.background = element_rect(fill = "white", color = "black"), 
    panel.spacing = unit(0,"lines"),
    strip.text.y = element_text(angle = 0)
  )
)


########### Pink
pink_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Pink Even"), 
             mapping = aes(x = cumPos, y = Fst, color = as_factor(chr)),
             alpha = 0.8, size = 0.7) +
  geom_point(data = filter(lrrc9_df, Species == "Pink Even"), mapping = aes(x = cumPos, y = Fst), 
             color = "#CC79A7", alpha = 1, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.006, 0.006)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.01)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02)) +
  ylab("Pink") +
  theme(plot.margin = margin(0.05,0.05,0,0.05, "in"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

########### SOCKEYE
sock_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Sockeye"), 
             mapping = aes(x = cumPos, y = Fst, color = as_factor(chr)),
             alpha = 0.8, size = 0.7) +
  geom_point(data = filter(lrrc9_df, Species == "Sockeye"), mapping = aes(x = cumPos, y = Fst), 
             color = "#CC79A7", alpha = 1, size = 0.8) +
  geom_point(data = filter(esr1_df, Species == "Sockeye"), mapping = aes(x = cumPos, y = Fst), 
             color = "#009E73", alpha = 1, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.006, 0.006)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02),
                     expand = expansion(mult = c(0.03, 0.01))) +
  ylab("Sockeye") +
  theme(plot.margin = margin(0,0.05,0,0.05, "in"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

########### CHUM
chum_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Chum"), 
             mapping = aes(x = cumPos, y = Fst, color = as_factor(chr)),
             alpha = 0.8, size = 0.7) +
  geom_point(data = filter(esrb_df, Species == "Chum"), mapping = aes(x = cumPos, y = Fst), 
             color = "#56B4E9", alpha = 1, size = 0.8) +
  geom_point(data = filter(lrrc9_df, Species == "Chum"), mapping = aes(x = cumPos, y = Fst), 
             color = "#CC79A7", alpha = 1, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.006, 0.006)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02),
                     expand = expansion(mult = c(0.03, 0.01))) +
  ylab("Chum") +
  theme(plot.margin = margin(0,0.05,0,0.05, "in"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) 

########### COHO
coho_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Coho"), 
             mapping = aes(x = cumPos, y = Fst, color = as_factor(chr)),
             alpha = 0.8, size = 0.7) +
  geom_point(data = filter(esrb_df, Species == "Coho"), mapping = aes(x = cumPos, y = Fst), 
             color = "#56B4E9", alpha = 1, size = 0.8) +
  geom_point(data = filter(esr1_df, Species == "Coho"), mapping = aes(x = cumPos, y = Fst), 
             color = "#009E73", alpha = 1, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.006, 0.006)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(0, 1.02),
                     expand = expansion(mult = c(0.03, 0.02))) +
  ylab("Coho") + xlab("Chromosome") +
  theme(plot.margin = margin(0,0.05,0.05,0.05, "in")) 


##################### COMBINE PLOTS ########################################

# add FST as separate label (it will be it's own plot)
y_lab <- ggplot() +
  annotate(geom = "text", x = 1, y = 1, angle = 90,
           label = expression(italic(F[ST])),
           size = mediumFS / .pt) + # need to convert annotation to pt 
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0),unit="in") # hopefully this removes whitespace around FST

######  ADD LABELS 
fourspp_cowplot <- plot_grid(pink_plot, sock_plot, chum_plot, coho_plot,
                             rel_heights = c(1,1,1,1.16),
                             ncol = 1,
                             labels = c('A','B','C', 'D'), 
                             label_fontfamily = "sans",
                             label_size = largeFS, label_colour = "black")

# combine FST label to other plots 
fourspp_cowplot_lab <- (y_lab - fourspp_cowplot) + # patchwork uses hyphen to allow for lefthand additions
  plot_layout(widths = c(1, 30))

png(paste0("./figures/fst/Figure1_wholegenome_trimlowFST_CBcolors_",format(Sys.Date(),"%Y%m%d"),".png"), 
    width = width2col, height = 0.7*width2col, res = 500, units = "in")
print(fourspp_cowplot_lab)
dev.off()

# i will try PDF
pdf(paste0("./figures/fst/Figure1_Wholegenome_trimlowFST_CBcolors_",format(Sys.Date(),"%Y%m%d"),".pdf"), 
    width = width2col, height = 0.7*width2col)
print(fourspp_cowplot_lab)
dev.off()

# this file size was ridiculous.... 
# svglite(str_c("./figures/fst/Figure1_wholegenome_",format(Sys.Date(),"%Y%m%d"),".svg"),
#         width= width2col, height= 0.7*width2col)
# print(fourspp_cowplot_lab)
# dev.off()
