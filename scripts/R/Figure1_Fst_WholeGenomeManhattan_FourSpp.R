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
pink_list <- as.list(Sys.glob("./results/fst/*pink-chum*minInd0.3*"))

# read in all data files that match wildcard
pink_df <- pink_list %>%
  set_names(nm = pink_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL, col.names = c("region", "chrName", "midPos", "Nsites", "Fst"))
  ) %>% 
  select(-c(region, Nsites)) %>%
  mutate(Species = "Pink Even")

########## SOCKEYE
sock_list <- as.list(Sys.glob("./results/fst/*euclide*minInd0.3*"))

# read in all data files that match wildcard
sock_df <- sock_list %>%
  set_names(nm = sock_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL, col.names = c("region", "chrName", "midPos", "Nsites", "Fst"))
  ) %>% 
  select(-c(region, Nsites)) %>%
  mutate(Species = "Sockeye")

########## CHUM
chum_list <- as.list(Sys.glob("./results/fst/*chumrun*minInd0.3*"))

# read in all data files that match wildcard
chum_df <- chum_list %>%
  set_names(nm = chum_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL, col.names = c("region", "chrName", "midPos", "Nsites", "Fst"))
  ) %>% 
  select(-c(region, Nsites)) %>%
  mutate(Species = "Chum")

########## COHO
coho_list <- as.list(Sys.glob("./results/fst/*coho-chum*minInd0.3*"))

# read in all data files that match wildcard
coho_df <- coho_list %>%
  set_names(nm = coho_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL, col.names = c("region", "chrName", "midPos", "Nsites", "Fst"))
  ) %>% 
  select(-c(region, Nsites)) %>%
  mutate(Species = "Coho")

rm(pink_list, sock_list, chum_list, coho_list)

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

rm(pink_df, sock_df, chum_df, coho_df) # remove to save space

# find max of each chrom
cumulate <- four_df %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(max_bp = max(midPos)) %>%
  mutate(final_bp = max_bp + 5) %>%
  mutate(cum_bp = cumsum(final_bp))

# create cumulative position column for all chromosomes
cum_df <- NULL
for(i in 1:length(unique(four_df$chr))){
  chrom <- cumulate$chr[i]
  # for chr == 1 only, bc not adding to any chrom.
  if(i == 1){cum_temp <- four_df %>%
    filter(chr == chrom) %>%
    mutate(cumPos = midPos)
  cum_df <- cum_temp
  }
  # for all other chromosomes, append
  if(i != 1){cumPos <- cumulate$cum_bp[i-1]
  cum_temp <- four_df %>%
    filter(chr == chrom) %>%
    mutate(cumPos = midPos + cumPos)
  cum_df <- rbind(cum_df,cum_temp)
  }
}

rm(cum_temp, chrom, i, cumPos) # remove unneeded files

# SINGLE OUT ESRB AND LRRC9 REGIONS
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

# Set the ggplot theme
theme_set(
  theme( 
    legend.position = "none",
    panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(angle = 0, size = 22),
    axis.title.x = element_text(size = 24, 
                                margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 26, angle = 90, 
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
             alpha = 0.8, size = 1.3) +
  geom_point(data = filter(lrrc9_df, Species == "Pink Even"), mapping = aes(x = cumPos, y = Fst), 
             color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02)) +
  ylab("Pink") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

########### SOCKEYE
sock_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Sockeye"), 
             mapping = aes(x = cumPos, y = Fst, color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  geom_point(data = filter(lrrc9_df, Species == "Sockeye"), mapping = aes(x = cumPos, y = Fst), 
             color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02),
                     expand = expansion(mult = c(0.02, 0.001))) +
  ylab("Sockeye") +
  theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

########### CHUM
chum_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Chum"), 
             mapping = aes(x = cumPos, y = Fst, color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  geom_point(data = filter(esrb_df, Species == "Chum"), mapping = aes(x = cumPos, y = Fst), 
             color = "dodgerblue2", alpha = 1, size = 1.5) +
  geom_point(data = filter(lrrc9_df, Species == "Chum"), mapping = aes(x = cumPos, y = Fst), 
             color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02),
                     expand = expansion(mult = c(0.02, 0.001))) +
  ylab("Chum") +
  theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) 

########### COHO
coho_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Coho"), 
             mapping = aes(x = cumPos, y = Fst, color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  geom_point(data = filter(esrb_df, Species == "Coho"), mapping = aes(x = cumPos, y = Fst), 
             color = "dodgerblue2", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(0, 1.02),
                     expand = expansion(mult = c(0.02, 0.001))) +
  ylab("Coho") + xlab("Chromosome") +
  theme(plot.margin = unit(c(0,0.1,0.1,0.1), "cm")) 


##################### COMBINE PLOTS ########################################

# add FST as separate label (it will be it's own plot)
y_lab <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = expression(italic(F[ST])), angle = 90, size = 16) + 
  coord_cartesian(clip = "off") +
  theme_void() 

######  ADD LABELS 
fourspp_cowplot <- plot_grid(pink_plot, sock_plot, chum_plot, coho_plot,
                             rel_heights = c(1,1,1,1.16),
                             ncol = 1,
                             labels = c('A','B','C', 'D'), 
                             label_fontfamily = "Arial",
                             label_size = 30, label_colour = "black")

# combine FST label to other plots 
fourspp_cowplot_lab <- (y_lab - fourspp_cowplot) + # patchwork uses hyphen to allow for lefthand additions
  plot_layout(widths = c(0.8, 20))

height = 14

jpeg(paste0("./figures/fst/fourspp_wholegenome_h",round(height, digits = 0),"_",format(Sys.Date(),"%Y%m%d"),".jpg"), 
     width = 20, height = height, res = 300, units = "in")
print(fourspp_cowplot_lab)
dev.off()


