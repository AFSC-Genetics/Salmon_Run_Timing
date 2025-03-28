# Create tree of comparisons between onohologs
# UPDATED 03/05/2025 with new trees with duplicated genes from all species

library(ape)
library(ggtree)
library(treeio)
library(tidyverse)

# Read the tree from Geneious
lrrc9_tree <- read.tree("./data/R/lrrc9_NJtree.newick")
plot(lrrc9_tree)

lrrc9LG <- data.frame(description = lrrc9_tree$tip.label,
                      tipname = gsub('_',' ',lrrc9_tree$tip.label),
                      lg = c("Outgroup","chr29","chr29","chr29","chr29",
                             "chr35","chr35","chr35","chr35")) %>%
  dplyr::mutate(FID = row_number(),
         tipname = gsub("'","",tipname),
         tipname = gsub(' LG.*','',tipname),
         tipname = gsub(' lrrc9','',tipname),
         tipname = gsub('L','',tipname))
  
lrrc9_tree[["node.label"]] <- round(as.numeric(lrrc9_tree[["node.label"]]), digits = 0)

lrrc9_homTree <- lrrc9_tree %>%
  as.treedata() %>%
  dplyr::mutate(FID = row_number()) %>%
  left_join(lrrc9LG, by = "FID") %>%
  treeio::rename_taxa(lrrc9LG, description, tipname) %>%
  ggtree(color = "black", size = 1) + 
  geom_tiplab(color = "black", size = 7) +
  geom_nodelab(hjust = 1, vjust=0, nudge_x = -0.004, nudge_y = 0.06, size = 6) +
  geom_treescale(x = 0.05, y = 1.2, linesize = 1, fontsize = 6) +
  theme(legend.position = "none",
        plot.margin = margin(1,1,1,1)) + 
  ggplot2::xlim(-0.02, 0.36)
lrrc9_homTree

# check nodes
lrrc9_homTree + geom_text(aes(label=node), color = "firebrick", size = 6)

lrrc9_homTree2 <- lrrc9_homTree %>%
  flip(9, 8) %>% # flip pink & chum lrrc9
  flip(4, 5) %>% # flip pink & coho lrrc9L
  flip(2, 13) %>% # flip chumlrrc9L w/ other lrrc9L
  flip(3, 14)     # flip socklrrc9L w/ other lrrc9L
lrrc9_homTree2

ggsave(paste0("./figures/tree/consensus/fourspp_lrrc9_consensus_tree_black_",format(Sys.Date(),"%Y%m%d"),".jpg"),
       lrrc9_homTree2, width = 8, height = 6)

####### ESRB ##################
# Read the tree from Geneious

esrb_tree <- read.tree("./data/R/esrb_NJtree.newick")
plot(esrb_tree)

esrb_tree[['tip.label']]

#### NEED TO MAKE CHANGES TO NEWICK FILE FOR THIS TO WORK PROPERLY
esrbLG <- data.frame(longname = esrb_tree$tip.label,
                     lg = c('Outgroup','chr35','chr35','chr35','chr35',
                            'chr29','chr29','chr29','chr29')) %>%
  dplyr::mutate(FID = row_number(),
         tipname = gsub('_',' ',longname),
         tipname = gsub("'","",tipname),
         tipname = gsub(' LG.*','',tipname),
         tipname = gsub(" esbr","",tipname), 
         tipname = gsub(' esrb','',tipname),
         tipname = gsub(' esr2b','',tipname),
         tipname = gsub('L','',tipname))

esrb_tree[["node.label"]] <- round(as.numeric(esrb_tree[["node.label"]]), digits = 0)

esrb_homTree <- esrb_tree %>%
  as.treedata() %>%
  dplyr::mutate(FID = row_number()) %>%
  left_join(esrbLG, by = "FID") %>%
  treeio::rename_taxa(esrbLG, longname, tipname) %>%
  ggtree(color = "black", size = 1) + 
  geom_tiplab(color = "black", size = 7) +
  geom_nodelab(hjust = 1, vjust=0, nudge_x = -0.004, nudge_y = 0.06, size = 6) +
  geom_treescale(x = 0.05, y = 1.2, linesize = 1, fontsize = 6) +
  theme(legend.position = "none") + 
  ggplot2::xlim(-0.01, 0.8)
esrb_homTree

esrb_homTree + geom_text(aes(label=node), color = "firebrick", size = 6)

esrb_homTree2 <- esrb_homTree %>%
  flip(8, 9) # flip pink & chum esrb
esrb_homTree2

ggsave(paste0("./figures/tree/consensus/fourspp_esrb_consensus_tree_black_",format(Sys.Date(),"%Y%m%d"),".jpg"),
       esrb_homTree2, width = 8, height = 6)

# END OF UPDATED SCRIPT

######### OLDER VERSION FROM JANUARY/FEB 2025 ########################

# Read the phylogenetic tree from Geneious
lrrc9_tree <- read.tree("./data/R/Chum_Sockeye_lrrc9like_lrrc9_PacificSalmon_Pike_Alignment_FocalSp_consensus_tree.newick")
plot(lrrc9_tree)

jpeg(paste0("./figures/tree/consensus/lrrc9_consensus_tree_",format(Sys.Date(), "%Y%m%d"),".jpg"), 
     width = 4, height = 6, res = 500, units = "in")
plot.phylo(
  lrrc9_tree,
  type = "phylogram",      # Tree layout
  align.tip.label = F,  # Align tip labels
  cex = 0.8,               # Font size for tip labels
  no.margin = TRUE         # Remove extra margins
)
dev.off()

# add data to say which chromosome it's on
lrrc9LG <- data.frame(description = lrrc9_tree$tip.label,
                      tipname = gsub("'","",lrrc9_tree$tip.label),
                      lg = c("Outgroup","chr29","chr29",
                             "chr35","chr35","chr35","chr35")) %>%
  dplyr::mutate(FID = row_number(),
         tipname = gsub("LG.*","",tipname),
         tipname = gsub("_"," ",tipname), tipname = gsub("  "," ",tipname),
         tipname = gsub("9like","9-like",tipname))

lrrc9_tree[["node.label"]] <- c("","100","100","100") # round values

lrrc9_homTree <- lrrc9_tree %>%
  as.treedata() %>%
  mutate(FID = row_number()) %>%
  left_join(lrrc9LG, by = "FID") %>%
  treeio::rename_taxa(lrrc9LG, description, tipname) %>%
  ggtree(color = "black", size = 1) + 
  geom_tiplab(color = "black", size = 6) +
  geom_nodelab(hjust = 1, vjust=0, nudge_x = -0.004, nudge_y = 0.06, size = 6) +
  geom_treescale(x = 0.05, y = 1.2, linesize = 1, fontsize = 6) +
  theme(legend.position = "none",
        plot.margin = margin(1,1,1,1)) + 
  ggplot2::xlim(-0.01, 0.36)
lrrc9_homTree



####### ESRB ##################
# Read the phylogenetic tree from Geneious
#esrb_tree2 <- read.tree("./data/R/esrb_tree.newick")
#  plot(esrb_tree)

esrb_tree <- read.tree("./data/R/esrb_align_esrbL_align_Alignment_FocalSp_consensus_tree.newick")
  plot(esrb_tree)
  
#### NEED TO MAKE CHANGES TO NEWICK FILE FOR THIS TO WORK PROPERLY
esrbLG <- data.frame(longname = esrb_tree$tip.label,
                     lg = c('Outgroup','chr35','chr35','chr35',
                            'chr29','chr29','chr29','chr29','chr29')) %>%
  mutate(FID = row_number(),
         tipname = gsub("LG.*","",longname),
         tipname = gsub("'","",tipname),
         #tipname = paste0(tipname,' (',lg,')'),
         tipname = gsub("_"," ",tipname), tipname = gsub("  "," ",tipname),
         tipname = gsub("esbr","esrb",tipname), tipname = gsub("bL","b-like",tipname),
         #tipname = gsub('\\s*\\(Outgroup\\)','',tipname),
         )

esrb_tree[["node.label"]] <- c("","100","100","93","100","100","96","96") # round values

esrb_homTree <- esrb_tree %>%
  as.treedata() %>%
  mutate(FID = row_number()) %>%
  left_join(esrbLG, by = "FID") %>%
  treeio::rename_taxa(esrbLG, longname, tipname) %>%
  ggtree(color = "black", size = 1) + 
  geom_tiplab(color = "black", size = 6) +
  geom_nodelab(hjust = 1, vjust=0, nudge_x = -0.004, nudge_y = 0.06, size = 6) +
  geom_treescale(x = 0.05, y = 1.2, linesize = 1, fontsize = 6) +
  theme(legend.position = "none") + 
  ggplot2::xlim(-0.01, 0.8)
esrb_homTree

ggsave(paste0("./figures/tree/consensus/fourspp_esrb_consensus_tree_black_",format(Sys.Date(),"%Y%m%d"),".jpg"),
       esrb_homTree, width = 8, height = 6)

##############################
unique(esrbLG$lg)
palette <- c("blue3","firebrick","gray")

esrb_homTree <- esrb_tree %>%
  as.treedata() %>%
  mutate(FID = row_number()) %>%
  left_join(esrbLG, by = "FID") %>%
  treeio::rename_taxa(esrbLG, longname, shortname_chr) %>%
  ggtree(aes(color = lg), size = 1) + 
  scale_color_manual(name = "Consensus Chum Chr", values = palette) +
  geom_tiplab(color = "black", size = 6) +
  geom_treescale(x = 0.05, y = 1.2, linesize = 1, fontsize = 7) +
  theme(legend.position = "none") + 
  ggplot2::xlim(0, 0.9)
esrb_homTree

ggsave(paste0("./figures/tree/consensus/esrb_consensus_tree_",format(Sys.Date(),"%Y%m%d"),".jpg"),
       esrb_homTree, width = 8, height = 10)

########### SCRATCH ############################################################
# COMBINE

jpeg(paste0("./figures/tree/consensus/both_consensus_trees_",format(Sys.Date(), "%Y%m%d"),".jpg"), 
     width = 10, height = 6, res = 500, units = "in")
# combine lrrc9 and esrb
plot_grid(plot.phylo(lrrc9_tree,type = "phylogram",align.tip.label = F,cex = 0.8,no.margin = TRUE), 
        NULL,
        plot.phylo(esrb_tree,type = "phylogram",align.tip.label = F,cex = 0.8,no.margin = TRUE),
        rel_widths = c(2,0.1,2), nrow = 1, ncol = 3,
        align = 'h',
        labels = c('B','','C'), 
        label_fontfamily = "helvetica",
        label_size = 18, label_colour = "black")
dev.off()

jpeg(paste0("./figures/tree/consensus/both_consensus_trees_",format(Sys.Date(), "%Y%m%d"),".jpg"), 
     width = 10, height = 6, res = 500, units = "in")
# combine lrrc9 and esrb
plot_grid(as.list(plot.phylo(lrrc9_tree,type = "phylogram",align.tip.label = F,cex = 0.8,no.margin = TRUE),
        plot.phylo(esrb_tree,type = "phylogram",align.tip.label = F,cex = 0.8,no.margin = TRUE)),
        rel_widths = c(2,2), nrow = 1, ncol = 3,
        align = 'h',
        labels = c('B','C'), #label_fontfamily = "Arial",
        label_fontfamily = "helvetica",
        label_size = 18, label_colour = "black"
)
dev.off()



homTree <- lrrc9_homTree2 + plot_spacer() + esrb_homTree +
  plot_layout(widths = c(1,0.4,1))
homTree







###### KIM CODE NEWICK FILE ########################

library(ape)
# Read the phylogenetic tree
#
tree <- read.tree("tree.newick")
plot(tree)
id <- read.csv("/home/kimberly.ledger/rockfish_mb/data/rkfish_ref_dbs/rockfish_reference_db_534_20250117.csv") %>%
  rename(accession_id = Name)
# Get species order
species_order_tissue <- as.data.frame(tree$tip.label) %>%
  rename(accession_id = "tree$tip.label") %>%
  left_join(id) %>%
  mutate(Organism = ifelse(accession_id == "'NC_060709 PCR Product'", "Sebastes ciliatus", Organism)) %>%
  mutate(Organism = ifelse(accession_id == "'NC_060711 PCR Product'", "Sebastes variabilis", Organism)) %>%
  mutate(Organism_abbr = sub("^(\\w)\\w*\\s(\\w+)", "\\1. \\2", Organism))
tree$tip.label <- species_order_tissue$Organism

plot.phylo(
  tree,
  type = "phylogram",      # Tree layout
  align.tip.label = F,  # Align tip labels
  cex = 0.8,               # Font size for tip labels
  no.margin = TRUE         # Remove extra margins
)
