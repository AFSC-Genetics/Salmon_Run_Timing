#https://nvelden.github.io/geneviewer/articles/geneviewer.html
#https://nvelden.github.io/geneviewer/articles/Examples.html


# NO LONGER INVERTING THE ZOOM IN


#library(devtools)
#devtools::install_github("nvelden/geneviewer")
library(geneviewer)
library(tidyverse)


wgd_genes <- read.csv("./data/R/wgd_geneviewer_data.csv", row.names = NULL) %>%
  filter(duplication != 7) # this is the middle ground value for the top plot

# REVERSE CHR35 SO ITS NO LONGER INVERTED
wgd_rev35 <- wgd_genes %>%
  mutate(start2 = start,
         end2 = end,
         start = ifelse(chr=='Oket35',-1*end2,start2),
         end = ifelse(chr=='Oket35',-1*start2,end2))

# PLOT WITH LABELS
GC_chart(wgd_rev35, 
         cluster = "chr",
         group = "duplication",
         height = 360) %>%
  GC_labels("Gene", cluster = 1,
            x = 0, y = 56, #dy = "-1.2em",dx = "0em",
            rotate = -50, #adjustLabels = TRUE, 
            fontSize = 14, fontFamily = "arial", fontStyle = "italic"
  ) %>%
  GC_labels("Gene", cluster = 2,
            x = 0, y = 10, #dy = "-1.2em", #dx = "2em",
            rotate = -50,
            fontSize = 14, fontFamily = "arial", fontStyle = "italic"
  ) %>%
  GC_genes(marker_size = "large", marker = "box", markerHeight = 30,
           customColors = c("#7570B3","#E7298A","#666666","#96cdcd","#87CEFF","#1C86EE")
  ) %>%
  GC_links(group = "duplication", 
           label = FALSE,
           use_group_colors = T) %>%
  #GC_clusterLabel(fontStyle = "italic", fontSize = 16, fontFamily = "arial") %>%
  GC_legend(F) %>%
  GC_scaleBar(cluster = 1, x = 0, y = 168,
              labelStyle = list(labelPosition = "right", fontWeight = "bold",
                                fontSize = 15, fontFamily = "arial"),
              textPadding = 4,
              scaleBarLineStyle = list(stroke = "black", strokeWidth = 1.5),
              scaleBarTickStyle = list(stroke = "black", strokeWidth = 1.5),
              title = "850kb", scaleBarUnit = 850000) %>%
  GC_scaleBar(cluster = 2, x = 0, y = 10,
              labelStyle = list(labelPosition = "right",
                                fontSize = 15, fontWeight = "bold",
                                fontFamily = "arial"),
              textPadding = 4,
              scaleBarLineStyle = list(stroke = "black", strokeWidth = 1.5),
              scaleBarTickStyle = list(stroke = "black", strokeWidth = 1.5),
              title = "800kb", scaleBarUnit = 800000)

# NO GENE LABELS
GC_chart(wgd_rev35, 
         cluster = "chr",
         group = "duplication",
         height = 360) %>%
  GC_genes(marker_size = "large", marker = "box", markerHeight = 30,
           customColors = c("#7570B3","#E7298A","#666666","#96cdcd","#87CEFF","#1C86EE")
  ) %>%
  GC_links(group = "duplication", 
           label = FALSE,
           use_group_colors = T) %>%
  #GC_clusterLabel(fontStyle = "italic", fontSize = 16, fontFamily = "arial") %>%
  GC_legend(F) %>%
  GC_scaleBar(cluster = 1, x = 0, y = 168,
              labelStyle = list(labelPosition = "right", fontWeight = "bold",
                                fontSize = 15, fontFamily = "arial"),
              textPadding = 4,
              scaleBarLineStyle = list(stroke = "black", strokeWidth = 1.5),
              scaleBarTickStyle = list(stroke = "black", strokeWidth = 1.5),
              title = "850kb", scaleBarUnit = 850000) %>%
  GC_scaleBar(cluster = 2, x = 0, y = 10,
              labelStyle = list(labelPosition = "right",
                                fontSize = 15, fontWeight = "bold",
                                fontFamily = "arial"),
              textPadding = 4,
              scaleBarLineStyle = list(stroke = "black", strokeWidth = 1.5),
              scaleBarTickStyle = list(stroke = "black", strokeWidth = 1.5),
              title = "800kb", scaleBarUnit = 800000)

# NO GENE LABELS
# ADD RED BORDER --> removed
GC_chart(wgd_rev35, 
         cluster = "chr",
         group = "duplication",
         #style = list(border = "3px solid red"),
         height = 360) %>%
  GC_genes(marker_size = "large", marker = "box", markerHeight = 30,
           customColors = c("#7570B3","#E7298A","#666666","#96cdcd","#87CEFF","#1C86EE")
  ) %>%
  GC_links(group = "duplication", 
           label = FALSE,
           use_group_colors = T) %>%
  #GC_clusterLabel(fontStyle = "italic", fontSize = 16, fontFamily = "arial") %>%
  GC_legend(F) %>%
  GC_scaleBar(cluster = 1, x = 0, y = 168,
              labelStyle = list(labelPosition = "right", fontWeight = "bold",
                                fontSize = 15, fontFamily = "arial"),
              textPadding = 4,
              scaleBarLineStyle = list(stroke = "black", strokeWidth = 1.5),
              scaleBarTickStyle = list(stroke = "black", strokeWidth = 1.5),
              title = " ", scaleBarUnit = 850000) %>%
  GC_scaleBar(cluster = 2, x = 0, y = 10,
              labelStyle = list(labelPosition = "right",
                                fontSize = 15, fontWeight = "bold",
                                fontFamily = "arial"),
              textPadding = 4,
              scaleBarLineStyle = list(stroke = "black", strokeWidth = 1.5),
              scaleBarTickStyle = list(stroke = "black", strokeWidth = 1.5),
              title = " ", scaleBarUnit = 800000) %>%
  GC_grid(margin = list(left = 5, right = 50))




# SIMPLE VERSION WITH ONLY LRRC9 AND ESRB FOR TOP PLOT
wgd_genes2 <- read.csv("./data/R/wgd_geneviewer_data.csv", row.names = NULL)

GC_chart(wgd_genes2, 
         cluster = "chr",
         group = "duplication",
         height = 200,
         style = list(border = "5px solid red")) %>%
  GC_align(id_column = "Gene", id = "FAM241A-like", align = "right"
  ) %>%
  GC_genes(marker_size = "large", marker = "box", markerHeight = 50,
           customColors = c("#7570B3","#E7298A","#666666","#96cdcd","#87CEFF","#1C86EE")
  ) %>%
  GC_links(group = "duplication", 
           value1 = c(1,2,3,4,5,6),
           value2 = c(1,2,3,4,5,6),
           linkWidth = 8,
           label = FALSE,
           use_group_colors = T) %>%
  GC_legend(F) %>%
  GC_grid(margin = list(top = 10, bottom = 0, left = 15, right = 15))

wgd_genes3 <- read.csv("./data/R/wgd_geneviewer_data.csv", row.names = NULL)[c(2,6,8,12:14),]

GC_chart(wgd_genes3, 
         cluster = "chr",
         group = "duplication",
         #style = list(border = "5px solid red"),
         height = 100) %>%
  GC_align(id_column = "Gene", id = "FAM241A-like", align = "right"
  ) %>%
  GC_genes(show = T,
           marker_size = "small", marker = "box",
           customColors = c("#E7298A","#1C86EE","white")
  ) %>%
  GC_links(group = "duplication", 
           value1 = c(2,6),
           value2 = c(2,6),
           linkWidth = 12,
           label = FALSE,
           use_group_colors = T) %>%
  GC_legend(F) %>%
  GC_grid(margin = list(top = 0, bottom = 0, left = 15, right = 15))

