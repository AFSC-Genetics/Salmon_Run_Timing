# setwd("../../Salmon_Run_Timing")

# PICK CREEK SOCKEYE ALIGNED TO CHUM FOR RUNTIMING STUDY AT LRRC9

# Pick Ck early and late run timings
#########################################################
# set up R
packages_needed <- c("dplyr", "tidyverse", "ggplot2", "readxl", "RColorBrewer", "ggpubr", "stats", "here")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

###########################################################################################################################

#to run locally
FEATURE_NAME <- "Runtime" #column name (in the metadata file) that contains the categorical variable you want to compare depth across (batch, region, sex, etc)
PREFIX <- "pick-chum" #the prefix for the lcWGS run
SPECIES <- "Sockeye"

DEPTHFILE <- paste0("./results/depth/",PREFIX,"_depths.csv")
#allmeta_df <- read.csv("./data/raw/fourspecies_runtiming_metadata.csv", header = T)
full_metadata  <- read.csv("./data/raw/pick_creek_metadata.csv", header = T)

features_df <- full_metadata %>%
  dplyr::rename(sampleID = sample) %>%
  dplyr::select(c("sampleID", FEATURE_NAME))

just_depths <- read.delim(DEPTHFILE, sep = "\t", header = F, row.names = NULL)
  colnames(just_depths) <- c("sampleID", "mean_depth")
  
just_depths$mean_depth <- as.numeric(just_depths$mean_depth)

depths_df <- left_join(features_df, just_depths, by = "sampleID")
  colnames(depths_df) <- c("sampleID", "feature", "mean_depth")
  
###########################################################################################################################
# plot depth distribution, colored by feature of interest

depths_plot <- ggplot(data = depths_df, aes(x = reorder(sampleID, -mean_depth), y = mean_depth, fill = feature)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  xlab("individual") + ylab("mean depth") +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = 1) +
  scale_y_continuous(breaks = c(0.5, 1, 2, 3, 4, 5)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        axis.text.x = element_blank(), axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank())

depths_plot

ggsave(paste0("./figures/depth/",PREFIX,"-",FEATURE_NAME,"_mean_depths.jpg"), 
       width = 8, height = 5, units = "in", dpi = 300)


###########################################################################################################################

# compare depth distributions (with a series of t-tests)
keeplist_df0.4 <- depths_df[depths_df$mean_depth >= 0.4,]
blocklist_df0.4 <- depths_df[depths_df$mean_depth < 0.4,]

write.table(blocklist_df0.4[,1], "./data/",PREFIX,"_blocklist_0.5x.txt",
            sep = "\t", quote = F, row.names = F, col.names = F)

keeplist_df0.4 %>%
  dplyr::group_by(feature) %>%
  dplyr::summarize(n=n(), avg_depth = mean(mean_depth)) 

