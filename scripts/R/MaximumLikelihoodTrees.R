# #### Maximum likelihood trees from Geneious Alignments
library(ape)
library(ggtree)
library(treeio)
library(phangorn)
library(tidyverse)

#### Estrogen Receptor Tree ####

#Read in alignment
ER_align_aa <- read.phyDat("./data/fasta/Esr1_aa_alignment_test.fasta", format = "fasta", type = "AA")

# Run a model test across common amino acid models
mt_results <- modelTest(ER_align_aa, model = c("WAG", "JTT", "LG", "Dayhoff", 
                                               "cpREV", "mtmam", "mtArt", "MtZoa", 
                                               "mtREV24", "VT","RtREV", "HIVw", 
                                               "HIVb", "FLU", "Blosum62", 
                                               "Dayhoff_DCMut", "JTT_DCMut") , multicore = TRUE)

best_model <- mt_results$Model[which.min(mt_results$AICc)]

# Calculate a distance matrix using the JTT model
dm <- dist.ml(ER_align_aa, model = "JTT")

# Generate a Neighbor-Joining starting tree
starting_tree <- NJ(dm)

# Set up the baseline likelihood object
# k = 4 for G(4), inv = 0.2 sets the starting estimate for invariant sites (+I)
fit_init <- pml(starting_tree, data = ER_align_aa, model = "JTT", k = 4, inv = 0.2)

# Optimize the model parameters and tree topology
fit_ml <- optim.pml(fit_init, 
                    model = "JTT", 
                    optNni = TRUE, 
                    optBf = FALSE,     # FALSE ensures standard JTT frequencies (No +F)
                    optGamma = TRUE,   # TRUE optimizes the Gamma shape parameter (+G)
                    optInv = TRUE)     # TRUE optimizes the proportion of invariant sites (+I)

print(fit_ml)

# Define the outgroup
my_outgroups <- c("ACF16007.1 estrogen receptor [Branchiostoma floridae]",
                  "NP_001191648.1 estrogen receptor [Aplysia californica]")

# Root the tree and put it back into the fit object
unrooted_tree <- fit_ml$tree
rooted_tree <- root(unrooted_tree, outgroup = my_outgroups, resolve.root = TRUE)
fit_ml$tree <- rooted_tree

# Run bootstrap analysis
bs_trees <- bootstrap.pml(fit_ml, bs = 100, optNni = TRUE)

# Plot the ML tree and overlay the bootstrap support values onto the nodes
plotBS(fit_ml$tree, bs_trees, type = "phylogram")

# Save tree with bootstrap values mapped to it
final_tree <- plotBS(fit_ml$tree, bs_trees, type = "none")
write.tree(final_tree, file = "ml_aa_tree.nwk")
#### format the ML tree
ER_tree <- final_tree

#ER_tree$node.label <- gsub(pattern="\\.0$",replacement="",ER_tree$node.label)

ggtree(ER_tree) + geom_tiplab(align=TRUE, family='mono',
                              linetype = "dotted", linesize = .7)

taxa_to_drop <- c("BAF45381.1 estrogen receptor [Magallana gigas]",
                  "XP_032811970.2 estrogen receptor isoform X1 [Petromyzon marinus]",
                  "esr1_PikeLG18_XP_010880455.1 esr1 [organism=Esox lucius] [GeneID=105017493] [isoform=X1]",
                  "NP_001191648.1 estrogen receptor [Aplysia californica]")

ER_tree <- drop.tip(ER_tree,taxa_to_drop)

#Rename the leaves
d <- read.csv("./data/raw/TipLabels.csv") %>%
  as_tibble() %>%
  rename('label' = 'Name_Seq' ,
         'newlabel' = 'NewTipLab' )

p <- ggtree(ER_tree) %<+% d

p2 <- p + geom_tiplab(aes(label=newlabel),align=TRUE, family='mono',
                      linetype = "dotted", linesize = .7)

#look at the prelim tree for formatting
p + geom_tiplab(aes(label=newlabel)) +geom_text(aes(label=node), color = "firebrick", size = 6)

#add metadata for node labels
node_metadata <- data.frame(
  node = c(28, 33, 47, 30),
  Event = c("1R/2R", "Ss4R", "Ss4R", "Ts3R")
)

p <- ggtree(ER_tree) %<+% d %<+% node_metadata

final_plot <- p + 
  geom_tiplab(aes(label=newlabel,
                  fontface = ifelse(newlabel %in% c("Chum esr2b*","Coho esr2b*",
                                                    "Atlantic Cod esr2a*","Atlantic Herring esr2a*",
                                                    "Coho esr1*", "Sockeye esr1*"), "bold", "plain"))) + 
  geom_nodepoint(aes(color = Event, shape = Event, size = Event), na.rm = TRUE) +
  
  # 1. Add bootstrap labels to the internal nodes
  # We subset using !isTip to only label internal nodes, and filter out empty labels
  geom_text2(aes(label = label, subset = !isTip & !is.na(as.numeric(label))), 
             hjust = 1.2, vjust = -0.5, size = 3, color = "darkgray") +
  
  # 2. Keep your special event labels (moved slightly right so they don't overlap)
  geom_text2(aes(label = Event), hjust = 1.2, vjust = -1.5, na.rm = TRUE, size = 3.5, fontface = "bold") +
  scale_color_manual(
    name = "Duplication Event",
    values = c("1R/2R" = "#E41A1C", "Ss4R" = "#377EB8", "Ts3R" = "#4DAF4C"),
    breaks = c("1R/2R", "Ts3R", "Ss4R"),
    labels = c("1R/2R, ~ 450MYA", "Ts3R, ~325MYA", "Ss4R, ~106MYA"),
    na.translate = FALSE
  ) +
  scale_shape_manual(
    name = "Duplication Event",
    values = c("1R/2R" = 17, "Ss4R" = 15, "Ts3R" = 19),
    breaks = c("1R/2R", "Ts3R", "Ss4R"),
    labels = c("1R/2R, ~ 450MYA", "Ts3R, ~325MYA", "Ss4R, ~106MYA"),
    na.translate = FALSE
  ) + 
  scale_size_manual(
    name = "Duplication Event",
    values = c("1R/2R" = 4, "Ss4R" = 4, "Ts3R" = 4),
    breaks = c("1R/2R", "Ts3R", "Ss4R"),
    labels = c("1R/2R, ~ 450MYA", "Ts3R, ~325MYA", "Ss4R, ~106MYA"),
    na.translate = FALSE
  ) +
  theme_tree2()+
  theme_tree() +  
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    
    # 1. Place the legend inside the plot window (relative x, y coordinates)
    # 0.15 moves it significantly to the left; 0.2 keeps it near the bottom
    legend.position = c(0.15, 0.65),
    
    # 2. Add a transparent background so it blends seamlessly into the whitespace
    legend.background = element_blank(),
    legend.box.background = element_blank()
  ) +
  # Draw the clean, dedicated phylogenetic scale bar
  geom_treescale(
    x = 0,               
    y = -1,              
    width = 0.5,         
    color = "black", 
    fontsize = 3
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35)))
print(final_plot)


#### Lrrc9 ML Tree ####
# #### Maximum likelihood tree from alignment

library(ape)
library(ggtree)
library(treeio)
library(phangorn)
library(tidyverse)

lrrc9_align <- read.phyDat("./data/fasta/lrrc9_Salmon_Alignment.fasta", format = "fasta", type = "DNA")

# Run a model test across common DNA models
mt_results <- modelTest(lrrc9_align, model = c("JC", "F81", "K80", "HKY", 
                                                "TrNe", "TrN", "TPM1", "K81", 
                                                "TPM1u", "TPM2", "TPM2u", "TPM3", 
                                                "TPM3u", "TIM1e", "TIM1", "TIM2e", 
                                                "TIM2", "TIM3e", "TIM3", "TVMe", 
                                                "TVM", "SYM", "GTR"), multicore = TRUE)

best_model <- mt_results$Model[which.min(mt_results$AICc)]

# Calculate a distance matrix using a basic model for the starting tree
dm <- dist.ml(lrrc9_align, model = "JC69")

# Generate a Neighbor-Joining starting tree
starting_tree <- NJ(dm)

# Set up the baseline likelihood object (k=4 sets up the 4 Gamma categories)
fit_init <- pml(starting_tree, data = lrrc9_align, k = 4)

# Optimize the model parameters and tree topology
fit_ml <- optim.pml(fit_init, 
                    model = "GTR", 
                    optNni = TRUE, 
                    optBf = TRUE,      # TRUE ensures unequal base frequencies (required for GTR)
                    optGamma = TRUE,   # TRUE optimizes the Gamma shape parameter (+G)
                    optInv = FALSE)    # FALSE prevents estimating invariant sites (+I)

print(fit_ml)

# Pike as the outgroup
my_outgroups <- c("Pike_lrrc9_LG15_NC_047583.1:16166890-16195893 lrrc9 [organism=Esox lucius] [GeneID=105015570] [chromosome=15]")

# Root the tree and put it back into the fit object
unrooted_tree <- fit_ml$tree
rooted_tree <- root(unrooted_tree, outgroup = my_outgroups, resolve.root = TRUE)
fit_ml$tree <- rooted_tree

# Bootstrap the ML tree
bs_trees <- bootstrap.pml(fit_ml, bs = 100, optNni = TRUE)

# Plot the ML tree and overlay the bootstrap support values onto the nodes
plotBS(fit_ml$tree, bs_trees, type = "phylogram")

# Save tree with bootstrap values mapped to it
final_tree <- plotBS(fit_ml$tree, bs_trees, type = "none")
write.tree(final_tree, file = "lrrc9_MLtree.nwk")

# Test for within and among ohnolog divergence within salmonids (drop)
# use permanova to deal with non-independence of branchlengths

#remove pike from the final tree
final_tree <-read.newick("./lrrc9_MLtree.nwk")
pike_tip <- final_tree$tip.label[grepl("pike", final_tree$tip.label, ignore.case = TRUE)]
pruned_tree <- ape::drop.tip(final_tree, pike_tip)
lrrc9Dists <- ape::cophenetic.phylo(pruned_tree)

Ohnolog_groupings <- data.frame(
  Taxon = rownames(lrrc9Dists),
  Group = ifelse(grepl("lrrc9L", rownames(lrrc9Dists)), "lrrc9L", "lrrc9")
)

set.seed(907)
permanova_lrrc9 <- vegan::adonis2(as.dist(lrrc9Dists) ~ Group, data = Ohnolog_groupings, permutations = 2000)

# r^2 = 0.984 genetic distance explained by gene family split
# p=0.029

#Test for heterogeneity in mutation rates (away from pike) 
lrrc9Dists_full <- ape::cophenetic.phylo(final_tree)
pike_name <- rownames(lrrc9Dists_full)[grepl("pike", rownames(lrrc9Dists_full), ignore.case = TRUE)]
pike_distances <- lrrc9Dists_full[pike_name, ]

rate_df <- data.frame(
  Taxon = names(pike_distances),
  Distance_to_Pike = as.numeric(pike_distances),
  Group = ifelse(grepl("lrrc9L", names(pike_distances)), "lrrc9L", "lrrc9")
)

rate_df <- rate_df[rate_df$Taxon != pike_name, ]

aggregate(Distance_to_Pike ~ Group, data = rate_df, FUN = mean) #mean distance to pike

lrrc9_wilcox <- wilcox.test(Distance_to_Pike ~ Group, data = rate_df)

#### format the ML tree for figure
lrrc9_tree <- final_tree

ggtree(lrrc9_tree) + geom_tiplab(align=TRUE, family='mono',
                              linetype = "dotted", linesize = .7)

#Rename the leaves
d <- read.csv("./data/raw/TipLabels_lrrc9.csv") %>%
  as_tibble() %>%
  rename('label' = 'Name_Seq' ,
         'newlabel' = 'NewTipLab' )

p <- ggtree(lrrc9_tree) %>% 
  rotate(17) %>% # put lrrc9 on top
  rotate(14) %>% # put coho on top 
  rotate(16) %>% # put sockeye under coho
  rotate(13) %>%
  rotate(12) %<+% d 

final_lrrc9_plot <- p + 
  geom_tiplab(aes(label=newlabel),size=6) + 
  #geom_nodelab(size = 3, hjust = -0.1, vjust = -0.3, color = "black") +
  # Draw the clean, dedicated phylogenetic scale bar
  geom_treescale(
    x = 0.15,               
    y = 1.15,              
    width = 0.05,         
    color = "black", 
    fontsize = 5,
    offset = 0.05
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35)))

png(filename = './figures/lrrc9_tree.png',width=5,height=5,units="in",res=300)
final_lrrc9_plot
dev.off()

#### esr2b ML Tree ####
# #### Maximum likelihood tree from alignment

library(ape)
library(ggtree)
library(treeio)
library(phangorn)
library(tidyverse)

esr2b_align <- read.phyDat("./data/fasta/esr2b_Salmon_Alignment.fasta", format = "fasta", type = "DNA")

# Run a model test across common DNA models
mt_results <- modelTest(esr2b_align, model = c("JC", "F81", "K80", "HKY", 
                                               "TrNe", "TrN", "TPM1", "K81", 
                                               "TPM1u", "TPM2", "TPM2u", "TPM3", 
                                               "TPM3u", "TIM1e", "TIM1", "TIM2e", 
                                               "TIM2", "TIM3e", "TIM3", "TVMe", 
                                               "TVM", "SYM", "GTR"), multicore = TRUE)

best_model <- mt_results$Model[which.min(mt_results$AICc)]

# Calculate a distance matrix using a basic model for the starting tree
dm <- dist.ml(esr2b_align, model = "JC69")

# Generate a Neighbor-Joining starting tree
starting_tree <- NJ(dm)

# Set up the baseline likelihood object (k=4 sets up the 4 Gamma categories)
fit_init <- pml(starting_tree, data = esr2b_align, k = 4)

# Optimize the model parameters and tree topology
fit_ml <- optim.pml(fit_init, 
                    model = "GTR", 
                    optNni = TRUE, 
                    optBf = TRUE,      # TRUE ensures unequal base frequencies (required for GTR)
                    optGamma = TRUE,   # TRUE optimizes the Gamma shape parameter (+G)
                    optInv = FALSE)    # FALSE prevents estimating invariant sites (+I)

print(fit_ml)

# Pike as the outgroup
my_outgroups <- c("Pike_esr2b_LG15_NC_047583.1:c16781261-16764829 esr2b [organism=Esox lucius] [GeneID=105015549] [chromosome=15]")

# Root the tree and put it back into the fit object
unrooted_tree <- fit_ml$tree
rooted_tree <- root(unrooted_tree, outgroup = my_outgroups, resolve.root = TRUE)
fit_ml$tree <- rooted_tree

# Bootstrap the ML tree
bs_trees <- bootstrap.pml(fit_ml, bs = 100, optNni = TRUE)

# Plot the ML tree and overlay the bootstrap support values onto the nodes
plotBS(fit_ml$tree, bs_trees, type = "phylogram")

# Save tree with bootstrap values mapped to it
final_tree <- plotBS(fit_ml$tree, bs_trees, type = "none")
write.tree(final_tree, file = "esr2b_MLtree.nwk")

# Test for within and among ohnolog divergence within salmonids (drop)
# use permanova to deal with non-independence of branchlengths

#remove pike from the final tree
final_tree <-read.newick("./esr2b_MLtree.nwk")
pike_tip <- final_tree$tip.label[grepl("pike", final_tree$tip.label, ignore.case = TRUE)]
pruned_tree <- ape::drop.tip(final_tree, pike_tip)

esr2bDists <- ape::cophenetic.phylo(final_tree)

Ohnolog_groupings <- data.frame(
  Taxon = rownames(esr2bDists),
  Group = ifelse(grepl("esrbL", rownames(esr2bDists)), "esr2bL", "esr2b") )%>% #note that coho is annotated esrb on both chromosomes. esr2bL ==14 and esr2b === 7
  mutate(Group = case_when(
                  grepl("Coho", Taxon) & grepl("LG14", Taxon) ~ "esr2bL",
                  TRUE ~ Group))%>%
  mutate(Group = case_when(
    grepl("Sockeye", Taxon) & grepl("LG18", Taxon) ~ "esr2b", #note that sockeye is annotated esrbL on both chromosomes. esr2bL ==12 and esr2b === 18
    TRUE ~ Group)) %>%
  mutate(Group = case_when(
    grepl("Pink", Taxon) & grepl("LG17", Taxon) ~ "esr2b", #note that pink is annotated esrbL on both chromosomes. esr2bL == 10 and esr2b === 17
    TRUE ~ Group)) %>%
  mutate(Group = case_when(
    grepl("Pink", Taxon) & grepl("LG10", Taxon) ~ "esr2bL", #note that pink is annotated esrbL on both chromosomes. esr2bL == 10 and esr2b === 17
    TRUE ~ Group))

set.seed(520)
permanova_esr2b <- vegan::adonis2(as.dist(esr2bDists) ~ Group, data = Ohnolog_groupings, permutations = 2000)

# r^2 = 0.93 of genetic distance is explained by ohnolog split
# p = 0.007

#Test for heterogeneity in mutation rates (away from pike) 
esr2bDists_full <- ape::cophenetic.phylo(final_tree)
pike_name <- rownames(esr2bDists_full)[grepl("pike", rownames(esr2bDists_full), ignore.case = TRUE)]
pike_distances <- esr2bDists_full[pike_name, ]

rate_df <- data.frame(
  Taxon = names(pike_distances),
  Distance_to_Pike = as.numeric(pike_distances),
  Group = ifelse(grepl("esrbL", names(pike_distances)), "esr2bL", "esr2b")) %>% #note that coho is annotated esrb on both chromosomes. esr2bL ==14 and esr2b === 7
  mutate(Group = case_when(
    grepl("Coho", Taxon) & grepl("LG14", Taxon) ~ "esr2bL",
    TRUE ~ Group))%>%
  mutate(Group = case_when(
    grepl("Sockeye", Taxon) & grepl("LG18", Taxon) ~ "esr2b", #note that sockeye is annotated esrbL on both chromosomes. esr2bL ==12 and esr2b === 18
    TRUE ~ Group)) %>%
  mutate(Group = case_when(
    grepl("Pink", Taxon) & grepl("LG17", Taxon) ~ "esr2b", #note that pink is annotated esrbL on both chromosomes. esr2bL == 10 and esr2b === 17
    TRUE ~ Group))%>%
  mutate(Group = case_when(
    grepl("Pink", Taxon) & grepl("LG10", Taxon) ~ "esr2bL", #note that pink is annotated esrbL on both chromosomes. esr2bL == 10 and esr2b === 17
    TRUE ~ Group))
  
rate_df <- rate_df[rate_df$Taxon != pike_name, ]

aggregate(Distance_to_Pike ~ Group, data = rate_df, FUN = mean) #mean distance to pike

esr2b_wilcox <- wilcox.test(Distance_to_Pike ~ Group, data = rate_df)

#p=0.02857
#W=0

#### format the ML tree
esr2b_tree <- final_tree

ggtree(esr2b_tree) + geom_tiplab(align=TRUE, family='mono',
                                 linetype = "dotted", linesize = .7)

#Rename the leaves
d <- read.csv("./data/raw/TipLabels_esr2b.csv") %>%
  as_tibble() %>%
  rename('label' = 'Name_Seq' ,
         'newlabel' = 'NewTipLab' )

p <- ggtree(esr2b_tree) %>% 
  rotate(13) %>%
  rotate(12) %>%
  rotate(16) %>%
  rotate(14) %>%
  rotate(17) %<+% d 

final_esr2b_plot <- p + 
  geom_tiplab(aes(label=newlabel),size=6) + 
  #geom_nodelab(size = 3, hjust = -0.1, vjust = -0.3, color = "black") +
  # Draw the clean, dedicated phylogenetic scale bar
  geom_treescale(
    x = 0.15,               
    y = 1.15,              
    width = 0.05,         
    color = "black", 
    fontsize = 5,
    offset = 0.05
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35)))

png(filename = './figures/esr2b_tree.png',width=5,height=5,units="in",res=300)
final_esr2b_plot
dev.off()
