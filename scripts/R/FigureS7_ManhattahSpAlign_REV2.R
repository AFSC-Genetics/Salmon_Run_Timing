# Species Specific Alignments
library(tidyverse)

# change to cb palette
geneColors <- c("#009E73", "#56B4E9", "#CC79A7")
wong_colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### Pink Salmon ####
#pink_list <- as.list(Sys.glob("./results/fst/SpeciesAlign/Pink/PhenotypeBased/*\\.fst\\.txt"))
pink_list <- as.list(Sys.glob("./results/fst/SpeciesSpecificAlign/Pink/PhenotypeBased/AukePink*.fst.txt"))

# read in all data files that match wildcard
pink_df <- pink_list %>%
  set_names(nm = pink_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL)
  ) %>%
  mutate(Species = "Pink Even")%>%
  select(region,chr,Nsites,Species)%>%
  rename('chr' = 1,
         'midPos' = 2,
         'Fst' = 3)

# # find max of each chrom
cumulatePink <- pink_df %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(max_bp = max(midPos), .groups = "drop") %>%
  mutate(bp_offset = lag(cumsum(max_bp + 5), default = 0))

pink_df <- pink_df %>%
  left_join(cumulatePink %>%
              dplyr::select(chr, bp_offset), by = "chr") %>%
  mutate(cumPos = midPos + bp_offset,
         Fst = as.numeric(Fst)) %>%
  dplyr::select(-bp_offset)

axis_setPink <- pink_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos)) %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))

# For Science, consistent font sizing
width2col <- 7.25  # 2 col, in inches

largeFS <- 10 # Panel Labels
mediumFS <- 8  # Axis titles, legend titles
smallFS <- 6  # Axis text, annotations, captions

# Set the ggplot theme - FROM FIG1
theme_set(
  theme( 
    legend.position = "none",
    panel.grid.major.y = element_blank(), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = smallFS, angle = 0),
    axis.title.x = element_text(size = mediumFS, 
                                margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = mediumFS, angle = 90, 
                                margin = margin(t = 0, r = 5, b = 0, l = 0)),
    panel.background = element_rect(fill = "white", color = "black"), 
    panel.spacing = unit(0,"lines"),
    strip.text.y = element_text(angle = 0)
  )
)

pink_plot <- ggplot() +
  geom_point(data = filter(pink_df, Species == "Pink Even"),
             mapping = aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setPink$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)),
                     label = axis_setPink$chr,
                     breaks = axis_setPink$center) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02)) +
  ylab("Pink") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')

PinkHighFst <- pink_df %>%
  filter(Fst > 0.18)

### lrrc9 NC_060182.1 72673514 - 72875058
pink_df %>%
  filter(chr == "NC_060182.1" & midPos >=72500000 & midPos <=73000000)%>%
  ggplot() + 
  geom_point(aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setPink$chr)))) +
  ylab("Pink") +
  geom_vline(xintercept=c(72650536,72875058),lty=2,color="#CC79A7")


lrrc9_pink <- pink_df %>%
  filter(chr == "NC_060182.1" & midPos >=72650536 & midPos <=72875058)

##### Pink Manhattans ####

q95 <- quantile(pink_df$Fst, 0.95, na.rm = TRUE)

pink_WG <- pink_df %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))%>%
  drop_na(Fst) %>%
  filter(Fst >= q95 | (Fst < q95 & runif(n()) < 0.2) ) %>%
  ungroup() %>%
  ggplot() + 
  geom_point(aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8, pch=20) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setPink$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_setPink$LG[c(T,F)], 
                     breaks = axis_setPink$center[c(T,F)]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.05)) +
  ylab("Pink") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text = element_text(angle = 0, size = smallFS),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  geom_point(data = filter(lrrc9_pink, Species == "Pink Even"), mapping = aes(x = cumPos, y = as.numeric(Fst)), 
             color = "#CC79A7", alpha = 1, size = 0.8, pch = 20)


#### Sockeye Salmon ####
sockeye_list <- as.list(Sys.glob("./results/fst/SpeciesAlign/Sockeye//*\\.fst\\.txt"))
sockeye_list <- as.list(Sys.glob("./results/fst/SpeciesSpecificAlign/Sockeye/euclide*.fst.txt"))

# read in all data files that match wildcard
sockeye_df <- sockeye_list %>%
  set_names(nm = sockeye_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL)
  ) %>%
  mutate(Species = "Sockeye")%>%
  select(region,chr,Nsites,Species)%>%
  rename('chr' = 1,
         'midPos' = 2,
         'Fst' = 3)

cumulateSockeye <- sockeye_df %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(max_bp = max(midPos), .groups = "drop") %>%
  mutate(bp_offset = lag(cumsum(max_bp + 5), default = 0))

sockeye_df <- sockeye_df %>%
  left_join(cumulateSockeye %>%
              dplyr::select(chr, bp_offset), by = "chr") %>%
  mutate(cumPos = midPos + bp_offset,
         Fst = as.numeric(Fst)) %>%
  dplyr::select(-bp_offset)

axis_setSockeye <- sockeye_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos))


lrrc9_sockeye <- sockeye_df %>% 
  filter(chr == "NC_042546.1" &  midPos >= 41162738 & midPos <= 41448689)

esr1_sockeye <- sockeye_df %>% 
  filter(chr == "NC_042547.1" &  midPos >= 7437083 & midPos <= 7631797) # esr1 ChrLG13: 7,621,613 - 7,645,274

sock_plot <- ggplot() +
  geom_point(data = filter(sockeye_df, Species == "Sockeye"),
             mapping = aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setSockeye$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)),
                     label = axis_setSockeye$chr,
                     breaks = axis_setSockeye$center) +
  scale_y_continuous(breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02),
                     expand = expansion(mult = c(0.02, 0.001))) +
  ylab("Sockeye") +
  theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')

SockeyeHighFst <- sockeye_df %>%
  mutate(Fst = as.numeric(Fst))%>%
  filter(Fst > 0.3)


## 12 NC_042546.1 41162738-41448689 lrrc9
ggplot() + 
  geom_point(data = filter(sockeye_df, Species == "Sockeye", chr == "NC_042546.1",
                           midPos >41000000, midPos<41725000), 
             mapping = aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.7, pch=20) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setSockeye$chr)))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_vline(xintercept=c(41162738,41448689),lty=2,color="#CC79A7")+
  ylab("Fst")

## 13 NC_042547.1 7437083 - 7631797 syne1a, esr1
# esr1 = 7,621,613 - 7,645,274 solidly in the highest Fst area
# syne1a 7,384,050-7,620,202
# armt1 7646409 - 7653175

esr1 <- c(7621613:7645274)
syne1a <- c(7384050:7620202)
armt1 <- c(7646409:7653175)

Sockeye_Chr13 <- sockeye_df %>%
  filter(Species == "Sockeye", chr == "NC_042547.1",midPos >6900000, midPos< 7850000) %>%
  mutate(gene = case_when(
    midPos %in% esr1 ~ 'esr1',
    midPos %in% syne1a ~ 'syne1a',
    midPos %in% armt1 ~ 'armt1',
    .default = 'No Gene'
  ))
#   
# ggplot() + 
#   geom_point(data =Chum_Chr13 %>% filter(gene=='No Gene'), aes(x = midPos, y = as.numeric(Fst), color = as_factor(gene)),
#              alpha = 0.3, size = 0.7) +
#   geom_point(data =Chum_Chr13 %>% filter(gene!='No Gene'), aes(x = midPos, y = as.numeric(Fst), color = as_factor(gene)),
#              alpha = 0.8, size = 0.8) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
#   ylab("Fst")+
#   scale_color_manual(
#     name = "Gene",        # Changes the title of the legend
#     values = wong_colors[c(5,2:4,6,7)] # Maps the string to your desired color
#   )+
#   theme(legend.position = "right")

## 19 NC_042553.1 31992649 No genes
ggplot() + 
  geom_point(data = filter(sockeye_df, Species == "Sockeye", chr == "NC_042553.1",midPos >31000000, midPos< 33500000), 
             mapping = aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setSockeye$chr)))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  #geom_vline(xintercept=c(31992649,31992649),lty=2,color='red')+
  ylab("Fst")

## 22 NC_042556.1 60913726 - 60988801 androgen receptor-like, eda2r, ophn1
ggplot() + 
  geom_point(data = filter(sockeye_df, Species == "Sockeye", chr == "NC_042556.1",midPos >60800000, midPos< 61075000), 
             mapping = aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setSockeye$chr)))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_vline(xintercept=c(60903181,60988802),lty=2,color="#CC79A7")+
  ylab("Fst")

##### Sockeye Manhattans ####
axis_setSockeye <- sockeye_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos)) %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))

q95 <- quantile(sockeye_df$Fst, 0.95, na.rm = TRUE)
  
sockeye_WG <- sockeye_df %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))%>%
  drop_na(Fst) %>%
  filter(Fst >= q95 | (Fst < q95 & runif(n()) < 0.2) ) %>%
  ungroup() %>%
  ggplot() + 
  geom_point(aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8, pch=20) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setSockeye$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = gsub(pattern = '^9$','9a',axis_setSockeye$LG[c(T,F)]) %>%gsub(pattern = '^29$','9b',.) , 
                     breaks = axis_setSockeye$center[c(T,F)]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.05)) +
  ylab("Sockeye") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text = element_text(angle = 0, size = smallFS),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  geom_point(data = filter(lrrc9_sockeye), mapping = aes(x = cumPos, y = as.numeric(Fst)), 
             color = "#CC79A7", alpha = 1, size = 0.8, pch=20)+
  geom_point(data = filter(esr1_sockeye), mapping = aes(x = cumPos, y = as.numeric(Fst)), 
             color = "#009E73", alpha = 1, size = 0.8, pch=20)


#### Coho Salmon ####
coho_list <- as.list(Sys.glob("./results/fst/SpeciesAlign/Coho/PhenotypeBased/*minDepthHalf*\\.fst\\.txt"))
coho_list <- as.list(Sys.glob("./results/fst/SpeciesSpecificAlign/Coho/PhenotypeBased/tocoho*minDepthHalf*.fst.txt"))

# read in all data files that match wildcard
coho_df <- coho_list %>%
  set_names(nm = coho_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL)
  ) %>%
  mutate(Species = "coho")%>%
  select(region,chr,Nsites,Species)%>%
  dplyr::rename('chr' = 1,
         'midPos' = 2,
         'Fst' = 3)

cumulateCoho <- coho_df %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(max_bp = max(midPos), .groups = "drop") %>%
  mutate(bp_offset = lag(cumsum(max_bp + 5), default = 0))

coho_df <- coho_df %>%
  left_join(cumulateCoho %>%
              dplyr::select(chr, bp_offset), by = "chr") %>%
  mutate(cumPos = midPos + bp_offset,
         Fst = as.numeric(Fst)) %>%
  dplyr::select(-bp_offset)

axis_setCoho <- coho_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos)) %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))


coho_plot <- ggplot() +
  geom_point(data = filter(coho_df, Species == "coho"),
             mapping = aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setCoho$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)),
                     label = axis_setCoho$chr,
                     breaks = axis_setCoho$center) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02)) +
  ylab("coho") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')


#3 peaks
CohoHighFst <- coho_df %>%
  mutate(Fst = as.numeric(Fst))%>%
  filter(Fst >= 0.5)

## 7 NC_034180.2 27117379 - 46013128 esr2b syne2b
ggplot() + 
  geom_point(data = filter(coho_df, Species == "coho", chr == "NC_034180.2",midPos >26000000, midPos<48000000), 
             mapping = aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setCoho$chr)))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_vline(xintercept=c(29455765,30509378),lty=2,color='#fe5a1d')+
  ylab("Fst")

## 17 NC_034190.2 CD9 antigen, ngfa (nerve growth factor)
ggplot() + 
  geom_point(data = filter(coho_df, Species == "coho", chr == "NC_034190.2",midPos >86000000, midPos< 88124368), 
             mapping = aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setCoho$chr)))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_vline(xintercept=c(86934443,87124368),lty=2,color='#fe5a1d')+
  ylab("Fst")

## 21 NC_034194.2 24744571 - 24902281 (syne1a), esr1, dcph1, zinc finger and BTB domain containing 2b, A kinase (PRKA) anchor protein 12b, methylenetetrahydrofolate dehydrogenase (NADP+ dependent) 1 like ## esr1 ChrLG21: 24,719,119 - 24,743,990 
#esr1 24,719,119 - 24,743,990
wong_colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
esr1 <- c(24719119:24743990)
dcph1 <- c(24745069:24754047)
LOC116356176 <- c(24759589:24760439)
zbtb2b <- c(24760519,24768737)
akap12b <- c(24770661:24803072)
mthfd1l <- c(24831113:24921764)

Coho_Chr7 <- coho_df %>%
  filter(Species == "coho", chr == "NC_034194.2",midPos >24300000, midPos< 25500000) %>%
  mutate(gene = case_when(
    midPos %in% esr1 ~ 'esr1',
    midPos %in% dcph1 ~ 'dcph1',
    midPos %in% LOC116356176 ~ 'LOC116356176',
    midPos %in% zbtb2b ~ 'zbtb2b',
    midPos %in% akap12b ~ 'akap12b',
    midPos %in% mthfd1l ~ 'mthfd1l',
    .default = 'No Gene'
  ))

ggplot() + 
  geom_point(data =Coho_Chr7 %>% filter(gene=='No Gene'), aes(x = midPos, y = as.numeric(Fst), color = as_factor(gene)),
             alpha = 0.3, size = 0.7) +
  geom_point(data =Coho_Chr7 %>% filter(gene!='No Gene'), aes(x = midPos, y = as.numeric(Fst), color = as_factor(gene)),
             alpha = 0.8, size = 0.8) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab("Fst")+
  scale_color_manual(
    name = "Gene",        # Changes the title of the legend
    values = wong_colors[c(5,2:4,6,7)] # Maps the string to your desired color
  )+
    theme(legend.position = "right")

# peak is from esr1 - mthfd1l
esr1[1]; mthfd1l[2]

esr1_coho <- coho_df %>%
  filter(chr == "NC_034194.2" & midPos >=esr1[1] & midPos <=mthfd1l[2]) 

esrb_coho <- coho_df %>%
  filter(chr == "NC_034180.2" & midPos >=29455765 & midPos <=30500057)



##### Coho Manhattans ####

q95 <- quantile(coho_df$Fst, 0.95, na.rm = TRUE)

coho_WG <- coho_df %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))%>%
  drop_na(Fst) %>%
  filter(Fst >= q95 | (Fst < q95 & runif(n()) < 0.2) ) %>%
  ungroup() %>%
  ggplot() + 
  geom_point(aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8, pch = 20) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setCoho$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_setCoho$LG[c(T,F)], 
                     breaks = axis_setCoho$center[c(T,F)]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.05)) +
  ylab("Coho") +
  xlab("Chromosome")+
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text = element_text(angle = 0, size = smallFS),
        axis.title.x = element_text(angle = 0, size = mediumFS),
        legend.position = 'none')+
  geom_point(data = filter(esr1_coho), mapping = aes(x = cumPos, y = as.numeric(Fst)), 
             color = "#009E73", alpha = 1, size = 0.8, pch = 20)+
  geom_point(data = filter(esrb_coho), mapping = aes(x = cumPos, y = as.numeric(Fst)), 
           color = "#56B4E9", alpha = 1, size = 0.8, pch = 20)



#### Chum Salmon ####
Chum_list <- as.list(Sys.glob("./results/fst/SpeciesAlign/Chum/*minDepthHalf*\\.fst\\.txt"))
Chum_list <- as.list(Sys.glob("./results/fst/SpeciesSpecificAlign/Chum/chumrun*fall-summer*minDepthHalf*.fst.txt"))

# read in all data files that match wildcard
Chum_df <- Chum_list %>%
  set_names(nm = Chum_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL)
  ) %>%
  mutate(Species = "Chum")%>%
  select(region,chr,Nsites,Species)%>%
  dplyr::rename('chr' = 1,
         'midPos' = 2,
         'Fst' = 3)

cumulateChum <- Chum_df %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(max_bp = max(midPos), .groups = "drop") %>%
  mutate(bp_offset = lag(cumsum(max_bp + 5), default = 0))

chum_df <- Chum_df %>%
  left_join(cumulateChum %>%
              dplyr::select(chr, bp_offset), by = "chr") %>%
  mutate(cumPos = midPos + bp_offset,
         Fst = as.numeric(Fst)) %>%
  dplyr::select(-bp_offset)

axis_setChum <- chum_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos)) %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))


Chum_plot <- ggplot() +
  geom_point(data = filter(chum_df, Species == "Chum"),
             mapping = aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setChum$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)),
                     label = axis_setChum$chr,
                     breaks = axis_setChum$center) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02)) +
  ylab("Chum") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')

ChumHighFst <- chum_df %>%
  filter(Fst > 0.4)

### LG15  NC_068435.1 18017912 - 19870544 tRNAs - bunch of genes
chum_df %>%
  filter(chr == "NC_068435.1" & midPos >=17500000 & midPos <=20500000)%>%
  ggplot() + 
  geom_point(aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setChum$chr)))) +
  ylab("Chum") +
  geom_vline(xintercept=c(18017912,19870544),lty=2,color="#CC79A7")

### LG23  NC_068443.1 10354488 - 10469110 kcnab1b
chum_df %>%
  filter(chr == "NC_068443.1" & midPos >=10000000 & midPos <=11000000)%>%
  ggplot() + 
  geom_point(aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setChum$chr)))) +
  ylab("Chum") +
  geom_vline(xintercept=c(10354488,10469110),lty=2,color="#CC79A7")

### LG28 NC_068448.1 850305 - 937353
chum_df %>%
  filter(chr == "NC_068448.1" & midPos >=750305 & midPos <=987353)%>%
  ggplot() + 
  geom_point(aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setChum$chr)))) +
  ylab("Chum") +
  geom_vline(xintercept=c(850305,937353),lty=2,color="#CC79A7")


### LG29  NC_068449.1 25440358 - 25545548 esrb
chum_df %>%
  filter(chr == "NC_068449.1" & midPos >=25000000 & midPos <=26545548)%>%
  ggplot() + 
  geom_point(aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setChum$chr)))) +
  ylab("Chum") +
  geom_vline(xintercept=c(25440358,25545548),lty=2,color="#CC79A7")

### LG35  NC_068455.1 28033205 - 28179022 lrrc9
chum_df %>%
  filter(chr == "NC_068455.1" & midPos >=27033205 & midPos <=28557254)%>%
  ggplot() + 
  geom_point(aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setChum$chr)))) +
  ylab("Chum") +
  geom_vline(xintercept=c(28033205,28206418),lty=2,color="#CC79A7")

esrb_chum <- chum_df %>%
  filter(chr == "NC_068449.1" & midPos >=25000000 & midPos <=26545548)

lrrc9_chum <- chum_df %>%
  filter(chr == "NC_068455.1" & midPos >=27033205 & midPos <=28557254)

er1_df <- chum_df %>%
  filter(chr == "NC_068428.1" & midPos >=23525841 & midPos <= 23552881)

Dummy_Leg <- chum_df %>% 
  filter((chr == "NC_068449.1" & midPos >=25000000 & midPos <=26545548) | (chr =="NC_068455.1" & midPos >=27033205 & midPos <=28557254) | (chr == "NC_068428.1" & midPos >=23525841 & midPos <= 23552881))%>%
  mutate(Color = ifelse(chr == "NC_068449.1" & midPos >=25000000 & midPos <=26545548,'esrb',
                        ifelse(chr == "NC_068455.1" & midPos >=27033205 & midPos <=28557254, 'lrrc9',
                               ifelse(chr == "NC_068428.1" & midPos >=23525841 & midPos <= 23552881,"esr1",NA))))%>%
  ggplot() + 
  geom_point(aes(x = cumPos, y = as.numeric(Fst), color = Color),
             alpha = 0.8, size = 0.8) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_setChum$LG[c(T,F)], 
                     breaks = axis_setChum$center[c(T,F)]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.05)) +
  ylab("Chum")+
  scale_color_manual(values = c("#009E73","#56B4E9","#CC79A7"),
                     labels = c('esr1','esr2b','lrrc9'),
                     name="Peak")+
  theme(legend.position = 'right')

p_leg <- cowplot::get_legend(Dummy_Leg) 


##### Chum Manhattans ####

q95 <- quantile(chum_df$Fst, 0.95, na.rm = TRUE)

Chum_WG <- chum_df %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))%>%
  drop_na(Fst) %>%
  filter(Fst >= q95 | (Fst < q95 & runif(n()) < 0.2) ) %>%
  ungroup() %>%
  ggplot() + 
  geom_point(aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 0.8, pch = 20) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_setChum$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_setChum$LG[c(T,F)], 
                     breaks = axis_setChum$center[c(T,F)]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.05)) +
  ylab("Chum") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text = element_text(angle = 0, size = smallFS),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  geom_point(data = lrrc9_chum, mapping = aes(x = cumPos, y = as.numeric(Fst)), 
             color = "#CC79A7", alpha = 1, size = 0.8, pch = 20)+
  geom_point(data = esrb_chum, mapping = aes(x = cumPos, y = as.numeric(Fst)), 
           color = "#56B4E9", alpha = 1, size = 0.8, pch = 20)


DummyY <- ggplot() +
  labs(y = expression(italic(F['ST']))) + # Set the desired y-axis label
  theme(axis.line.x =  element_blank(),      # Rm the x-axis line
    axis.title.x =  element_blank(),     # Rm the x-axis title
    axis.text.x =  element_blank(),      # Rm the x-axis text (tick labels)
    #axis.title.y = element_blank(),    # Keep the y-axis title
    axis.text.y = element_blank(),     # Remove the y-axis text
    axis.ticks.y = element_blank(),    # Remove the y-axis ticks
    legend.position = "none",          # Remove the legend
    panel.background = element_blank(),# Remove the panel background
    panel.border = element_blank(),    # Remove the panel border
    panel.grid.major = element_blank(),# Remove major grid lines
    panel.grid.minor = element_blank(),# Remove minor grid lines
    plot.background = element_blank()  # Remove the plot background
  )


#Combine plots
pSpcAlgn <- cowplot::plot_grid(DummyY,
                   cowplot::plot_grid(pink_WG,
                                      sockeye_WG,
                                      Chum_WG,
                                      coho_WG,
                                      nrow=4,
                                      labels = c("A", "B", "C", "D"),
                                      label_fontfamily = "sans",
                                      label_size = largeFS,
                                      label_colour = "black",
                                      rel_widths = c(98,98,98,100)),
                   #p_leg,
                   nrow=1,
                   ncol=2, #ncol=3, # w/o legend
                   rel_widths = c(3,98)) #,100)) # w/o legend

png("./figures/FigureS7_SpeciesSpecificAlignment_trimLowFst.png",
    height=width2col,width=1.4*width2col,res=300,units='in')
pSpcAlgn
dev.off()


pdf("./figures/FigureS7_SpeciesSpecificAlignment_trimLowFst.pdf",
    height=width2col,width=1.4*width2col)
pSpcAlgn
dev.off()



