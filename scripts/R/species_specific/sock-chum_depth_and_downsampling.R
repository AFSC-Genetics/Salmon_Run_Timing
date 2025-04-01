# SOCKEYE ALIGNED TO CHUM FOR RUNTIMING STUDY AT LRRC9

# Whitefish Ck males (Late Stream) compared to Teal Ck (Early, from Euclide et al., 2024)
# they have different depths of coverage, but also different samples sizes
# after looking at this paper [https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14286],
# decided that the effective sample size is not drastically different between two pops
# this is when all jacks are removed from whitefish data

# code originally written by Laura Timm
# NH edited for project

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
PREFIX <- "sock-chum" #the prefix for the lcWGS run
SPECIES <- "Sockeye"

DEPTHFILE <- paste0("./results/depth/",PREFIX,"_depths.csv")
METADATAFILE <- "./data/raw/fourspecies_runtiming_metadata.csv"

###########################################################################################################################

# NH edited Laura Timm's depth function
depths_t_test <- function(df, colnm, f1, f2) {
  first_factor_depths <- df[which(df[[colnm]] == f1), which(colnames(df) == "mean_depth")]
  second_factor_depths <- df[which(df[[colnm]] == f2), which(colnames(df) == "mean_depth")]
  if (length(first_factor_depths) < 5) {
    print(paste0("WARNING: ", f1, " is represented by fewer than five individuals. Consider removing these from the depths analysis and adding ", f1, " to the downsample list (factors_to_downsample) manually."))
  }
  if (length(second_factor_depths) < 5) {
    print(paste0("WARNING: ", f2, " is represented by fewer than five individuals. Consider removing these from the depths analysis and adding ", f2, " to the downsample list (factors_to_downsample) manually."))
  }
  f1_depths = as.numeric(first_factor_depths)
  f2_depths = as.numeric(second_factor_depths)
  t_out <- t.test(x = f1_depths, y = f2_depths, var.equal = FALSE)
  outlist <- list(pval = t_out$p.value, f1_depths = f1_depths, f2_depths = f2_depths)
  return(outlist)
}

###########################################################################################################################

full_metadata <- read.csv(METADATAFILE) %>%
  filter(Species == SPECIES)
head(full_metadata)

features_df <- full_metadata[c("sampleID", FEATURE_NAME)]
  head(features_df)

just_depths <- read.delim(DEPTHFILE, sep = "\t", header = F, row.names = NULL)
  colnames(just_depths) <- c("sampleID", "mean_depth")
  head(just_depths)

just_depths$mean_depth <- as.numeric(just_depths$mean_depth)

depths_df <- left_join(features_df, just_depths, by = "sampleID")
  colnames(depths_df) <- c("sampleID", "feature", "mean_depth")
  head(depths_df)

# remove whitefish samples that weren't ever sequenced but are in the metadata
depths_df <- filter(depths_df, !is.na(mean_depth))

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
keeplist_df0.5 <- depths_df[depths_df$mean_depth >= 0.5,]
blocklist_df0.5 <- depths_df[depths_df$mean_depth < 0.5,]

write.table(blocklist_df0.5[,1], "./data/",PREFIX,"_blocklist_0.5x.txt",
            sep = "\t", quote = F, row.names = F, col.names = F)

###########################################################################
######### DID NOT RUN THE FOLLOWING BC DECIDED NOT TO DOWNSAMPLE ##########
###########################################################################
keeplist_df <- keeplist_df0.5
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
downsampled_factors

# No downsampling because of effective population size of more Whitefish samples than Euclide samples