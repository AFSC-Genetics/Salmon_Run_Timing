# PINK ALIGNED TO CHUM FOR RUNTIMING STUDY AT LRRC9

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
PREFIX <- "pink-chum" #the prefix for the lcWGS run
SPECIES = "Pink"

DEPTHFILE <- paste0("./data/depth/",PREFIX,"_depths.csv")
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
        axis.text.x = element_blank(), axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank())

depths_plot

ggsave(paste0("./figures/depth/",PREFIX,"-",FEATURE_NAME,"_mean_depths.jpg"), width = 8, height = 5, units = "in", dpi = 300)

###########################################################################################################################

# Remove samples below 0.5x, which there were none
keeplist_df0.5 <- depths_df[depths_df$mean_depth >= 0.5,]
blocklist_df0.5 <- depths_df[depths_df$mean_depth < 0.5,]

write.table(blocklist_df0.5[,1], paste0("./data/R/",PREFIX,"blocklist_0.5x.txt"),
            sep = "\t", quote = F, row.names = F, col.names = F)

####### CHECK IF DOWNSAMPLING SHOULD BE CONSIDERED ########################

#BASED ON CHOSEN CUTOFF, WHAT KEEPLIST ARE YOU USING
keeplist_df <- keeplist_df0.5
keeplist_df <- keeplist_df1

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
  downsampled_factors
}


######## MANUALLY DOWNSAMPLE THE TOP 13 SAMPLES (SINCE THERE WAS A JUMP IN DEPTH THERE)

depths_order <- depths_df[order(depths_df$mean_depth, decreasing = T),]
depths_downsample <- depths_order[1:13,]

theoretical_depths <- depths_downsample[1:13,] %>% 
  mutate(mean_depth = 2.6)

depths_check_downsample <- rbind(theoretical_depths, depths_order[14:nrow(depths_order),])

for (i in 1:length(factors_list)) {
  for (j in 2:length(factors_list)) {
    if (i < j) {
      first_factor <- as.character(factors_list[i])
      second_factor <- as.character(factors_list[j])
      results <- depths_t_test(depths_check_downsample, "feature", first_factor, second_factor)
      if (results$pval < 0.05) {
          print("Need to select more individuals")  
        }else{  
          print("This will correct for difference in depth")
        }
     }
  }
}

# first 13 get new depths
depths_new <- depths_downsample %>%
  mutate(seed_depth = 42 + (2.6 / mean_depth)) %>%
  mutate(array_input = paste0(row_number(),":",sampleID,":",round(seed_depth,digits = 3)))

write.table(depths_new[,"array_input"], file = paste0("./data/",PREFIX,"_13indivs_downsampleARRAY_input.txt"),
            row.names = F, col.names = F, quote = F)

##########################################################################################################
# RERAN MEAN_DEPTH.PY SCRIPT TO GET NEW DEPTHS AND PLOT AFTER DOWNSAMPLING

post_down_depth <- read.csv(paste0("./data/depth/",PREFIX,"_downsample_depths.csv"), header = F,
                            col.names = c("sampleID", "new_depth")) 

post_down_all_depths <- full_join(depths_order, post_down_depth, by = "sampleID")

post_down_all_depths <- post_down_all_depths %>%
  mutate(new_depth = ifelse(is.na(new_depth), mean_depth, new_depth))

###### Re-Plot

#plot new depth values, if needed
downdepths_plot <- ggplot(data = post_down_all_depths, aes(x = reorder(sampleID, -new_depth), y = new_depth, fill = feature)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  xlab("individual") +
  ylab("mean depth") +
  ggtitle(paste0("Pink-Chum: Downsampled")) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = 1) +
  scale_y_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5),
                     expand = expansion(0.01,0.01)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


ggsave(paste0("./figures/depth/",PREFIX, "-", FEATURE_NAME, "_mean_downsampled_depths.jpg"), 
       plot = downdepths_plot, width = 8, height = 5, units = "in", dpi = 300)
