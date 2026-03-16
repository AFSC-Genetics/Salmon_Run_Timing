# plot WHOLE GENOME pink-odd and sockeye (Teal & Whitefish)

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

###### Read in each per-chrom Fst File ########

########## PINK
pink_list <- as.list(Sys.glob("../2024_pinkOdd/results/fst/*pink-odd*minInd0.3*fst.txt"))

# read in all data files that match wildcard
pink_df <- pink_list %>%
  set_names(nm = pink_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL, col.names = c("region", "chrName", "midPos", "Nsites", "Fst"))
  ) %>% 
  select(-c(region, Nsites)) %>%
  mutate(Species = "Pink - Odd")

########## SOCKEYE - Whitefish-Teal
sock_list <- as.list(Sys.glob("../2024_sockeye/results/fst/sock-chum*early-late*minInd0.3*fst.txt"))

# read in all data files that match wildcard
sock_df <- sock_list %>%
  set_names(nm = sock_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL, col.names = c("region", "chrName", "midPos", "Nsites", "Fst"))
  ) %>% 
  select(-c(region, Nsites)) %>%
  mutate(Species = "Sockeye - Streams")

########## SOCKEYE - Whitefish-Teal
pick_list <- as.list(Sys.glob("../2024_sockeye/results/fst/pick-chum*early-late*minInd0.3*fst.txt"))

# read in all data files that match wildcard
pick_df <- pick_list %>%
  set_names(nm = pick_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL, col.names = c("region", "chrName", "midPos", "Nsites", "Fst"))
  ) %>% 
  select(-c(region, Nsites)) %>%
  mutate(Species = "Sockeye - Pick Creek")

rm(pink_list, sock_list, pick_list)

############## COMBINE THREE FST OUTPUTS INTO ONE DATAFRAME #######################
# bind rows together
comb_df <- bind_rows(pink_df, sock_df, pick_df) %>%
  mutate(midPos = as.numeric(midPos/10^6),
         Fst = as.numeric(Fst) %>% ifelse(. < 0,0,.))

#add chromosome numbers, remove chrName
comb_df <- left_join(comb_df, chrom_df, by = "chrName") %>%
  select(chr, midPos, Fst, Species)

rm(pink_df, sock_df, pick_df) # remove to save space

# find max of each chrom
cumulate <- comb_df %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(max_bp = max(midPos)) %>%
  mutate(final_bp = max_bp + 5) %>%
  mutate(cum_bp = cumsum(final_bp))

# create cumulative position column for all chromosomes
cum_df <- NULL
for(i in 1:length(unique(comb_df$chr))){
  chrom <- cumulate$chr[i]
  # for chr == 1 only, bc not adding to any chrom.
  if(i == 1){cum_temp <- comb_df %>%
    filter(chr == chrom) %>%
    mutate(cumPos = midPos)
  cum_df <- cum_temp
  }
  # for all other chromosomes, append
  if(i != 1){cumPos <- cumulate$cum_bp[i-1]
  cum_temp <- comb_df %>%
    filter(chr == chrom) %>%
    mutate(cumPos = midPos + cumPos)
  cum_df <- rbind(cum_df,cum_temp)
  }
}

rm(cum_temp, chrom, i, cumPos) # remove unneeded files

lrrc9_df <- cum_df %>%
  filter(chr == 35,
         midPos >= 27.85,
         midPos <= 28.20)

esr1_df <- cum_df %>%
  filter(chr == 8,
         midPos >= 23.141735,
         midPos <= 23.906114)


############# PREPARE FOR PLOTTING #########################################

axis_set <- cum_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos)) %>%
  filter(chr %% 2 == 1)

# Set the ggplot theme
theme_set(
  theme( 
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(angle = 0, size = 22),
    axis.title.x = element_text(size = 24, 
                                margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 24, angle = 90, 
                                margin = margin(t = 0, r = 5, b = 0, l = 0)),
    panel.background = element_rect(fill = "white", color = "black"), 
    panel.spacing = unit(0,"lines"),
    strip.text.y = element_text(angle = 0)
  )
)

########### Pink
pink_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Pink - Odd"), 
             mapping = aes(x = cumPos, y = Fst, color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  geom_point(data = filter(lrrc9_df, Species == "Pink - Odd"), mapping = aes(x = cumPos, y = Fst), 
             color = "brown1", alpha = 1, size = 1.6) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02),
                     expand = expansion(mult = c(0.02, 0.001))) +
  ylab("Pink \n (Odd)") +
  theme(plot.margin = unit(c(0.6,0.1,0,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) 

########### Sockeye - Pick Creek
sock_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Sockeye - Streams"), 
             mapping = aes(x = cumPos, y = Fst, color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  geom_point(data = filter(lrrc9_df, Species == "Sockeye - Streams"), mapping = aes(x = cumPos, y = Fst), 
             color = "brown1", alpha = 1, size = 1.6) +
  geom_point(data = filter(esr1_df, Species == "Sockeye - Streams"), mapping = aes(x = cumPos, y = Fst), 
             color = "green3", alpha = 1, size = 1.6) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02),
                     expand = expansion(mult = c(0.02, 0.001))) +
  ylab("Sockeye \n (Whitefish-Teal)") +
  theme(plot.margin = unit(c(0.6,0.1,0,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) 

########### SOCKEYE - Whitefish v Teal
pick_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Sockeye - Pick Creek"), 
             mapping = aes(x = cumPos, y = Fst, color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  # geom_point(data = filter(lrrc9_df, Species == "Sockeye - Pick Creek"), mapping = aes(x = cumPos, y = Fst), 
  #            color = "brown1", alpha = 1, size = 1.6) +
  geom_point(data = filter(esr1_df, Species == "Sockeye - Pick Creek"), mapping = aes(x = cumPos, y = Fst), 
             color = "green3", alpha = 1, size = 1.6) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(0, 1.02),
                     expand = expansion(mult = c(0.02, 0.001))) +
  ylab("Sockeye \n (Pick Creek)") +
  xlab("Chromosome") +
  theme(plot.margin = unit(c(0.4,0.1,0.1,0.1), "cm"))

##################### COMBINE PLOTS ########################################

# add FST as separate label (it will be it's own plot)
y_lab <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = expression(italic(F[ST])), angle = 90, size = 16) + 
  coord_cartesian(clip = "off") +
  theme_void() 

###### ADD LABELS #####
threespp_cowplot <- plot_grid(pink_plot, sock_plot, pick_plot,
                            rel_heights = c(1, 1, 1.15),
                            ncol = 1,
                            labels = c('A','B','C'), 
                            label_fontfamily = "ArialMT",
                            label_x = 0, label_y = 1.06,
                            label_size = 30, label_colour = "black")

# combine FST label to other plots 
threespp_cowplot_lab <- (y_lab - threespp_cowplot) + # patchwork uses hyphen to allow for lefthand additions
  plot_layout(widths = c(0.8, 20))

height = 11

jpeg(paste0("./figures/fst/supplfig_threespp_wholegenome_LABEL2_h",round(height, digits = 0),"_",format(Sys.Date(),"%Y%m%d"),".jpg"), 
     width = 20, height = height, res = 300, units = "in")
print(threespp_cowplot_lab)
dev.off()


# test <- (y_lab - pick_plot) + plot_layout(widths = c(0.8, 20))
# 
# jpeg(paste0("./figures/fst/supplfig_test_h",round(height, digits = 0),"_",format(Sys.Date(),"%Y%m%d"),".jpg"),
#      width = 20, height = 5, res = 300, units = "in")
# print(test)
# dev.off()


