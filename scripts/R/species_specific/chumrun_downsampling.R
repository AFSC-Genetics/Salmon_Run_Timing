# Chum Run Timing Depths 

# depth/downsampling code originally written by Laura Timm
# Natasha Howe edited for project

#########################################################
# set up R
packages_needed <- c("dplyr", "tidyverse", "ggplot2", "readxl", "ggpubr", "stats", "here")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
} 
  rm(packages_needed, i)
  
###########################################################################################################################

#to run locally
FEATURE_NAME <- "Runtime" #column name (in the metadata file) that contains the categorical variable you want to compare depth across (batch, region, sex, etc)
PREFIX <- "chumrun" #the prefix for the lcWGS run
SPECIES <- "Chum" # for filtering the metadata table

setwd(here())
DEPTHFILE <- paste0("./data/depth/",PREFIX,"_depths.csv")
METADATAFILE <- "./data/raw/yukon_chum_metadata.csv"
METADATAFILE <- "./data/raw/fourspecies_runtiming_metadata.csv"
###########################################################################################################################

###########################################################################################################################

## ADD THE MATURITY DATA TO FILTER OUT JACKS
full_metadata <- read.csv(METADATAFILE) %>%
  filter(Species == SPECIES)
head(full_metadata)

features_df <- full_metadata[c("sampleID", FEATURE_NAME)]
head(features_df)

just_depths <- read.csv(DEPTHFILE, header = F, row.names = NULL)
colnames(just_depths) <- c("sampleID", "mean_depth")
head(just_depths)

just_depths$mean_depth <- as.numeric(just_depths$mean_depth)

depths_df <- left_join(just_depths, features_df, by = "sampleID")
colnames(depths_df) <- c("sampleID", "mean_depth", "feature")
head(depths_df)

###########################################################################################################################
# plot depth distribution, colored by feature of interest

depths_plot <- ggplot(data = depths_df, aes(x = reorder(sampleID, -mean_depth), y = mean_depth, fill = feature)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  xlab("individual") + ylab("mean depth") +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = 1) +
  scale_y_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5),
                     expand = expansion(0.01,0.01)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

depths_plot

ggsave(paste0(BASEDIR, "./figures/depth/",PREFIX, "-", FEATURE_NAME, "_mean_depths.jpg"), width = 8, height = 5, units = "in", dpi = 300)


###########################################################################################################################

# compare depth distributions (with a series of t-tests)

# This is our general standard for this low of depths 
keeplist_df0.5 <- depths_df[depths_df$mean_depth >= 0.5,]
blocklist_df0.5 <- depths_df[depths_df$mean_depth < 0.5,]

# However, because there is only two below 0.5 and they're not far off, we keep those.
keeplist_df0.4 <- depths_df[depths_df$mean_depth >= 0.39,] # similar including the one that rounded up to 0.4
blocklist_df0.4 <- depths_df[depths_df$mean_depth < 0.39,] # similar including the one that rounded up to 0.4

# check mean depth after 0.4x cutoff (none removed)
mean((filter(keeplist_df0.4, feature == "Early"))$mean_depth)
mean((filter(keeplist_df0.4, feature == "Late"))$mean_depth)

write.table(blocklist_df0.4[,1], "./sedna_files/inputs/blocklist_0.4x.txt",
            sep = "\t", quote = F, row.names = F, col.names = F)

###########################################################################
####### CHECK IF DOWNSAMPLING SHOULD BE CONSIDERED ########################

#BASED ON CHOSEN CUTOFF, WHAT KEEPLIST ARE YOU USING
keeplist_df <- keeplist_df0.4

factors_list <- unique(keeplist_df[["feature"]])

## get a list of factors (batch1, batch2, etc) that require downsampling
downsampled_factors <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(downsampled_factors) <- c("lo_factor", "hi_factor", "targeted_mean", "targeted_sd")

for (i in 1:length(factors_list)) {
  for (j in 2:length(factors_list)) {
    if (i < j) {
      first_factor <- as.character(factors_list[i])
      second_factor <- as.character(factors_list[j])
      results <- depths_t_test(keeplist_df, "feature", first_factor, second_factor)
      if (results$pval < 0.05) {
        if (mean(results$f1_depths) < mean(results$f2_depths)) {
          lo_mean <- mean(results$f1_depths)
          lo_sd <- sd(results$f1_depths)
          target_factor <- first_factor
          downsample_factor <- second_factor
        } else {
          lo_mean <- mean(results$f2_depths)
          lo_sd <- sd(results$f2_depths)
          target_factor <- second_factor
          downsample_factor <- first_factor
        }
        downsampled_factors <- rbind(downsampled_factors, 
                                     list("lo_factor" = as.character(target_factor),
                                          "hi_factor" = as.character(downsample_factor), 
                                          "targeted_mean" = as.double(lo_mean),
                                          "targeted_sd" = as.double(lo_sd)))
      }
    }
  }
}
#did it work?
if(nrow(downsampled_factors) == 0){
  print("Welch's t-test showed no significant difference in depths across early and late run timing groups.")
}else{
  print("Welch's t-test showed significant difference in depths across early and late run timing groups.")
}

