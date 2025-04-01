# Create tree of comparisons between onhologs

library(ape)
library(ggtree)
library(treeio)
library(tidyverse)

####### LRRC9 ##################################################################

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

####### ESRB ####################################################################
# Read the tree from Geneious

esrb_tree <- read.tree("./data/R/esrb_NJtree.newick")
plot(esrb_tree)

esrb_tree[['tip.label']]

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

