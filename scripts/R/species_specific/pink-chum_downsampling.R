# PINK ALIGNED TO CHUM FOR RUNTIMING STUDY AT LRRC9

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
PREFIX <- "pink-chum" #the prefix for the lcWGS run

BASEDIR <- here()
DEPTHFILE <- paste0("./results/depth/",PREFIX,"_depths.csv")
METADATAFILE <- "./data/raw/pink_even_collection_062024.csv"

###########################################################################################################################

# define function

# NH EDITED THIS TO WORK FOR MY DATASET
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
full_metadata <- read.csv(METADATAFILE)
head(full_metadata)

features_df <- full_metadata[c("ABLG", FEATURE_NAME)]
head(features_df)

just_depths <- read.csv(DEPTHFILE, header = F, row.names = NULL)
colnames(just_depths) <- c("ABLG", "mean_depth")
head(just_depths)

just_depths <- just_depths %>%
  mutate(ABLG = gsub("ABLG", "", ABLG))

just_depths$ABLG <- as.numeric(just_depths$ABLG)
just_depths$mean_depth <- as.numeric(just_depths$mean_depth)

depths_df <- left_join(just_depths, features_df, by = "ABLG")
colnames(depths_df) <- c("ABLG", "mean_depth", "feature")
head(depths_df)

###########################################################################################################################
# plot depth distribution, colored by feature of interest

depths_plot <- ggplot(data = depths_df, aes(x = reorder(ABLG, -mean_depth), y = mean_depth, fill = feature)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  xlab("individual") +
  ylab("mean depth") +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = 1) +
  scale_y_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5),
                     expand = expansion(0.01,0.01)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

depths_plot

ggsave(paste0(BASEDIR, "./figures/depth/",PREFIX, "-", FEATURE_NAME, "_mean_depths.jpg"), width = 8, height = 5, units = "in", dpi = 300)


###########################################################################################################################

# compare depth distributions (with a series of t-tests)

out_df <- depths_df %>%
  mutate(sampleID = paste0("ABLG",ABLG)) %>%
  select(sampleID, feature, mean_depth)
head(out_df)

keeplist_df1 <- out_df[out_df$mean_depth >= 1,]
blocklist_df1 <- out_df[out_df$mean_depth < 1,]

keeplist_df0.5 <- out_df[out_df$mean_depth >= 0.5,]
blocklist_df0.5 <- out_df[out_df$mean_depth < 0.5,]

write.table(blocklist_df0.5[,1], "./sedna_files/inputs/blocklist_0.5x.txt",
            sep = "\t", quote = F, row.names = F, col.names = F)

write.table(blocklist_df1[,1], "./sedna_files/inputs/blocklist_1x.txt",
            sep = "\t", quote = F, row.names = F, col.names = F)

# for FST input
write.table(keeplist_df0.5[,c(1:2)], "./sedna_files/inputs/early_late_keeplist0.5x_fst.txt",
            sep = "\t", quote = F, row.names = F, col.names = F)

###########################################################################
####### CHECK IF DOWNSAMPLING SHOULD BE CONSIDERED ########################

#BASED ON CHOSEN CUTOFF, WHAT KEEPLIST ARE YOU USING
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
  mutate(array_input = paste0(row_number(),":ABLG",ABLG,":",round(seed_depth,digits = 3)))

write.table(depths_new[,"array_input"], file = paste0("./sedna_files/inputs/",PREFIX,"_13indivs_downsampleARRAY_input.txt"),
            row.names = F, col.names = F, quote = F)

##########################################################################################################
# RERAN MEAN_DEPTH.PY SCRIPT TO GET NEW DEPTHS AND PLOT AFTER DOWNSAMPLING

post_down_depth <- read.csv(paste0("./results/depth/",PREFIX,"_downsample_depths.csv"), header = F,
                            col.names = c("ABLG", "new_depth")) 

post_down_depth <-post_down_depth %>%
  mutate(ABLG = gsub("ABLG", "", ABLG))

post_down_depth$ABLG <- as.numeric(post_down_depth$ABLG)

post_down_all_depths <- full_join(depths_order, post_down_depth, by = "ABLG")

post_down_all_depths <- post_down_all_depths %>%
  mutate(new_depth = ifelse(is.na(new_depth), mean_depth, new_depth))

# MAY NEED TO EDIT CODE MORE TO WORK FOR BELOW                    !!!!!!!!!!!!!!!!!!!!!!!!!

#plot new depth values, if needed
downdepths_plot <- ggplot(data = post_down_all_depths, aes(x = reorder(ABLG, -new_depth), y = new_depth, fill = feature)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  xlab("individual") +
  ylab("mean depth") +
  ggtitle(paste0("Pink-Chum: Downsampled 13 individuals to 2.6x")) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = 1) +
  scale_y_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5),
                     expand = expansion(0.01,0.01)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

downdepths_plot

ggsave(paste0(BASEDIR, "./figures/depth/",PREFIX, "-", FEATURE_NAME, "_mean_downsampled_depths.jpg"), width = 8, height = 5, units = "in", dpi = 300)


################ DID NOT USE THIS BELOW ########################################
################################################################################

# write out subsetBAMS_array_input

if (exists("downsample_df")) {
  outfile <- file(paste0(here(), PREFIX, "-", FEATURE_NAME, "_downsample-bamsARRAY_inputTEST.txt"), "wb")
  iterator <- 1
  for (sample in 1:nrow(downsample_df)) {
    if (downsample_df[sample, "target_prop", 1] != 1) {
      outline <- paste0(as.character(iterator), ":", downsample_df[sample, "sample"], ":", as.character(round(as.numeric(downsample_df[sample, "target_prop"]), 3)))
      write(outline, file = outfile, append = TRUE)
      iterator <- iterator + 1
    }
  }
  
  close(outfile)
  
  #the file that was written is uploaded to sedna and used as the input array file for lcWGSpipeline_opt-downsample.py.
  #Be sure -a flag in opt-downsample matches the name of the file written here
  #unless you change it, the filename is "PREFIX-FEATURENAME_downsample-bamsARRAY_input.txt":
  print(paste0(PREFIX, "-", FEATURE_NAME, "_downsample-bamsARRAY_input.txt"))
}
