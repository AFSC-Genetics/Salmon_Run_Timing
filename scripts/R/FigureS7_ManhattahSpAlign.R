# Species Specific Alignments
library(tidyverse)


#### Pink Salmon ####
pink_list <- as.list(Sys.glob("./results/fst/SpeciesAlign/Pink/PhenotypeBased/*\\.fst\\.txt"))

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

# find max of each chrom
cumulate <- pink_df %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(max_bp = max(midPos)) %>%
  mutate(final_bp = max_bp + 5) %>%
  mutate(cum_bp = cumsum(final_bp))

# create cumulative position column for all chromosomes
cum_df <- NULL
for(i in 1:length(unique(pink_df$chr))){
  chrom <- cumulate$chr[i]
  # for chr == 1 only, bc not adding to any chrom.
  if(i == 1){cum_temp <- pink_df %>%
    filter(chr == chrom) %>%
    mutate(cumPos = midPos)
  cum_df <- cum_temp
  }
  # for all other chromosomes, append
  if(i != 1){cumPos <- cumulate$cum_bp[i-1]
  cum_temp <- pink_df %>%
    filter(chr == chrom) %>%
    mutate(cumPos = midPos + cumPos)
  cum_df <- rbind(cum_df,cum_temp)
  }
}

axis_set <- cum_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos))

# Set the ggplot theme
theme_set(
  theme( 
    legend.position = "none",
    panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(angle = 0, size = 22),
    axis.title.x = element_text(size = 18, 
                                margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 24, angle = 90, 
                                margin = margin(t = 0, r = 5, b = 0, l = 0)),
    panel.background = element_rect(fill = "white", color = "black"), 
    panel.spacing = unit(0,"lines"),
    strip.text.y = element_text(angle = 0)
  )
)

pink_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Pink Even"), 
             mapping = aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  #geom_point(data = filter(lrrc9_df, Species == "Pink Even"), mapping = aes(x = cumPos, y = Fst), 
  #           color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02)) +
  ylab("Pink") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')

PinkHighFst <- cum_df %>%
  filter(Fst > 0.18)

### lrrc9 NC_060182.1 72673514 - 72875058
cum_df %>%
  filter(chr == "NC_060182.1" & midPos >=72500000 & midPos <=73000000)%>%
ggplot() + 
  geom_point(aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  #geom_point(data = filter(lrrc9_df, Species == "Pink Even"), mapping = aes(x = cumPos, y = Fst), 
  #           color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  ylab("Pink") +
  geom_vline(xintercept=c(72650536,72875058),lty=2,color='#fe5a1d')

lrrc9_df <- cum_df %>%
  filter(chr == "NC_060182.1" & midPos >=72650536 & midPos <=72875058)

##### Pink Manhattans ####
axis_set <- cum_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos)) %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))

pink_WG <- cum_df %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))%>%
  ggplot() + 
  geom_point(aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$LG[c(T,F)], 
                     breaks = axis_set$center[c(T,F)]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.05)) +
  ylab("Pink") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text = element_text(angle = 0, size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  geom_point(data = filter(lrrc9_df, Species == "Pink Even"), mapping = aes(x = cumPos, y = as.numeric(Fst)), 
             color = "#fe5a1d", alpha = 1, size = 1.5)



#### Sockeye Salmon ####

sockeye_list <- as.list(Sys.glob("./results/fst/SpeciesAlign/Sockeye//*\\.fst\\.txt"))

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

# find max of each chrom
cumulate <- sockeye_df %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(max_bp = max(midPos)) %>%
  mutate(final_bp = max_bp + 5) %>%
  mutate(cum_bp = cumsum(final_bp))

# create cumulative position column for all chromosomes
cum_df <- NULL
for(i in 1:length(unique(sockeye_df$chr))){
  chrom <- cumulate$chr[i]
  # for chr == 1 only, bc not adding to any chrom.
  if(i == 1){cum_temp <- sockeye_df %>%
    filter(chr == chrom) %>%
    mutate(cumPos = midPos)
  cum_df <- cum_temp
  }
  # for all other chromosomes, append
  if(i != 1){cumPos <- cumulate$cum_bp[i-1]
  cum_temp <- sockeye_df %>%
    filter(chr == chrom) %>%
    mutate(cumPos = midPos + cumPos)
  cum_df <- rbind(cum_df,cum_temp)
  }
}

axis_set <- cum_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos))


# Set the ggplot theme
# Set the ggplot theme
theme_set(
  theme( 
    legend.position = "none",
    panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(angle = 0, size = 22),
    axis.title.x = element_text(size = 18, 
                                margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 24, angle = 90, 
                                margin = margin(t = 0, r = 5, b = 0, l = 0)),
    panel.background = element_rect(fill = "white", color = "black"), 
    panel.spacing = unit(0,"lines"),
    strip.text.y = element_text(angle = 0)
  )
)

lrrc9_df <- cum_df %>% 
  filter(chr == "NC_042546.1" &  midPos >= 41162738 & midPos <= 41448689)

esr1_df <- cum_df %>% 
  filter(chr == "NC_042547.1" &  midPos >= 7437083 & midPos <= 7631797) # esr1 ChrLG13: 7,621,613 - 7,645,274

sock_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Sockeye"), 
             mapping = aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  #geom_point(data = filter(lrrc9_df, Species == "Sockeye"), mapping = aes(x = cummidPos, y = Fst), 
  #           color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02),
                     expand = expansion(mult = c(0.02, 0.001))) +
  ylab("Sockeye") +
  theme(plot.margin = unit(c(0,0.1,0,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')

SockeyeHighFst <- cum_df %>%
  mutate(Fst = as.numeric(Fst))%>%
  filter(Fst > 0.3)

## 12 NC_042546.1 41162738-41448689 lrrc9
ggplot() + 
  geom_point(data = filter(cum_df, Species == "Sockeye", chr == "NC_042546.1",midPos >41000000, midPos<41725000), 
             mapping = aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_vline(xintercept=c(41162738,41448689),lty=2,color='#fe5a1d')+
  ylab("Fst")

## 13 NC_042547.1 7437083 - 7631797 syne1a, esr1
# esr1 = 7,621,613 - 7,645,274 solidly in the highest Fst area
# syne1a 7,384,050-7,620,202
# armt1 7646409 - 7653175

esr1 <- c(7621613:7645274)
syne1a <- c(7384050:7620202)
armt1 <- c(7646409:7653175)

Chum_Chr13 <- cum_df %>%
  filter(Species == "Sockeye", chr == "NC_042547.1",midPos >6900000, midPos< 7850000) %>%
  mutate(gene = case_when(
    midPos %in% esr1 ~ 'esr1',
    midPos %in% syne1a ~ 'syne1a',
    midPos %in% armt1 ~ 'armt1',
    .default = 'No Gene'
  )) 
  

  ggplot() + 
    geom_point(data =Chum_Chr13 %>% filter(gene=='No Gene'), aes(x = midPos, y = as.numeric(Fst), color = as_factor(gene)),
               alpha = 0.3, size = 1) +
    geom_point(data =Chum_Chr13 %>% filter(gene!='No Gene'), aes(x = midPos, y = as.numeric(Fst), color = as_factor(gene)),
               alpha = 0.8, size = 1.3) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    ylab("Fst")+
    scale_color_manual(
      name = "Gene",        # Changes the title of the legend
      values = wong_colors[c(5,2:4,6,7)] # Maps the string to your desired color
    )+
    theme(legend.position = "right")

## 19 NC_042553.1 31992649 No genes
ggplot() + 
  geom_point(data = filter(cum_df, Species == "Sockeye", chr == "NC_042553.1",midPos >31000000, midPos< 33500000), 
             mapping = aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  #geom_vline(xintercept=c(31992649,31992649),lty=2,color='red')+
  ylab("Fst")

## 22 NC_042556.1 60913726 - 60988801 androgen receptor-like, eda2r, ophn1
ggplot() + 
  geom_point(data = filter(cum_df, Species == "Sockeye", chr == "NC_042556.1",midPos >60800000, midPos< 61075000), 
             mapping = aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_vline(xintercept=c(60903181,60988802),lty=2,color='#fe5a1d')+
  ylab("Fst")

##### Sockeye Manhattans ####
axis_set <- cum_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos)) %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))

sockeye_WG <- cum_df %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))%>%
  ggplot() + 
  geom_point(aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = gsub(pattern = '^9$','9a',axis_set$LG[c(T,F)]) %>%gsub(pattern = '^29$','9b',.) , 
                     breaks = axis_set$center[c(T,F)]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.05)) +
  ylab("Sockeye") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text = element_text(angle = 0, size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.Position = 'none')+
  geom_point(data = filter(lrrc9_df), mapping = aes(x = cumPos, y = as.numeric(Fst)), 
             color = "#fe5a1d", alpha = 1, size = 1.5)+
  geom_point(data = filter(esr1_df), mapping = aes(x = cumPos, y = as.numeric(Fst)), 
             color = "#4CBB17", alpha = 1, size = 1.5)


#### Coho Salmon ####
coho_list <- as.list(Sys.glob("./results/fst/SpeciesAlign/Coho/PhenotypeBased/*minDepthHalf*\\.fst\\.txt"))

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

# find max of each chrom
cumulate <- coho_df %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(max_bp = max(midPos)) %>%
  mutate(final_bp = max_bp + 5) %>%
  mutate(cum_bp = cumsum(final_bp))

# create cumulative position column for all chromosomes
cum_df <- NULL
for(i in 1:length(unique(coho_df$chr))){
  chrom <- cumulate$chr[i]
  # for chr == 1 only, bc not adding to any chrom.
  if(i == 1){cum_temp <- coho_df %>%
    filter(chr == chrom) %>%
    mutate(cumPos = midPos)
  cum_df <- cum_temp
  }
  # for all other chromosomes, append
  if(i != 1){cumPos <- cumulate$cum_bp[i-1]
  cum_temp <- coho_df %>%
    filter(chr == chrom) %>%
    mutate(cumPos = midPos + cumPos)
  cum_df <- rbind(cum_df,cum_temp)
  }
}

axis_set <- cum_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos)) 

# Set the ggplot theme
theme_set(
  theme( 
    legend.position = "none",
    panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(angle = 0, size = 22),
    axis.title.x = element_text(size = 24, 
                                margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 26, angle = 90, 
                                margin = margin(t = 0, r = 5, b = 0, l = 0)),
    panel.background = element_rect(fill = "white", color = "black"), 
    panel.spacing = unit(0,"lines"),
    strip.text.y = element_text(angle = 0)
  )
)

coho_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "coho"), 
             mapping = aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  #geom_point(data = filter(lrrc9_df, Species == "coho Even"), mapping = aes(x = cumPos, y = Fst), 
  #           color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
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
CohoHighFst <- cum_df %>%
  mutate(Fst = as.numeric(Fst))%>%
  filter(Fst >= 0.5)

## 7 NC_034180.2 27117379 - 46013128 esr2b syne2b
ggplot() + 
  geom_point(data = filter(cum_df, Species == "coho", chr == "NC_034180.2",midPos >26000000, midPos<48000000), 
             mapping = aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_vline(xintercept=c(29455765,30509378),lty=2,color='#fe5a1d')+
  ylab("Fst")

## 17 NC_034190.2 CD9 antigen, ngfa (nerve growth factor)
ggplot() + 
  geom_point(data = filter(cum_df, Species == "coho", chr == "NC_034190.2",midPos >86000000, midPos< 88124368), 
             mapping = aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
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

Coho_Chr7 <- cum_df %>%
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
             alpha = 0.3, size = 1) +
  geom_point(data =Coho_Chr7 %>% filter(gene!='No Gene'), aes(x = midPos, y = as.numeric(Fst), color = as_factor(gene)),
             alpha = 0.8, size = 1.3) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab("Fst")+
  scale_color_manual(
    name = "Gene",        # Changes the title of the legend
    values = wong_colors[c(5,2:4,6,7)] # Maps the string to your desired color
  )+
    theme(legend.position = "right")

esrb_df <- cum_df %>%
  filter(chr == "NC_034180.2" & midPos >=29455765 & midPos <=30500057)

esr1_df <- cum_df %>%
  filter(chr == "NC_034194.2" & midPos >=24719119 & midPos <=24743990) #24,719,119..24,743,990

##### Coho Manhattans ####
axis_set <- cum_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos)) %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))

coho_WG <- cum_df %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))%>%
  ggplot() + 
  geom_point(aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$LG[c(T,F)], 
                     breaks = axis_set$center[c(T,F)]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.05)) +
  ylab("Coho") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text = element_text(angle = 0, size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  geom_point(data = filter(esr1_df), mapping = aes(x = cumPos, y = as.numeric(Fst)), 
             color = "#4CBB17", alpha = 1, size = 1.5)+
  geom_point(data = filter(esrb_df), mapping = aes(x = cumPos, y = as.numeric(Fst)), 
           color = "dodgerblue2", alpha = 1, size = 1.5)



#### Chum Salmon ####
Chum_list <- as.list(Sys.glob("./results/fst/SpeciesAlign/Chum/*minDepthHalf*\\.fst\\.txt"))

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

# find max of each chrom
cumulate <- Chum_df %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(max_bp = max(midPos)) %>%
  mutate(final_bp = max_bp + 5) %>%
  mutate(cum_bp = cumsum(final_bp))

# create cumulative position column for all chromosomes
cum_df <- NULL
for(i in 1:length(unique(Chum_df$chr))){
  chrom <- cumulate$chr[i]
  # for chr == 1 only, bc not adding to any chrom.
  if(i == 1){cum_temp <- Chum_df %>%
    filter(chr == chrom) %>%
    mutate(cumPos = midPos)
  cum_df <- cum_temp
  }
  # for all other chromosomes, append
  if(i != 1){cumPos <- cumulate$cum_bp[i-1]
  cum_temp <- Chum_df %>%
    filter(chr == chrom) %>%
    mutate(cumPos = midPos + cumPos)
  cum_df <- rbind(cum_df,cum_temp)
  }
}

axis_set <- cum_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos)) 

# Set the ggplot theme
theme_set(
  theme( 
    legend.position = "none",
    panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(angle = 0, size = 22),
    axis.title.x = element_text(size = 24, 
                                margin = margin(t = 5, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 26, angle = 90, 
                                margin = margin(t = 0, r = 5, b = 0, l = 0)),
    panel.background = element_rect(fill = "white", color = "black"), 
    panel.spacing = unit(0,"lines"),
    strip.text.y = element_text(angle = 0)
  )
)

Chum_plot <- ggplot() + 
  geom_point(data = filter(cum_df, Species == "Chum"), 
             mapping = aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  #geom_point(data = filter(lrrc9_df, Species == "Chum Even"), mapping = aes(x = cumPos, y = Fst), 
  #           color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$chr, 
                     breaks = axis_set$center) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.02)) +
  ylab("Chum") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')

ChumHighFst <- cum_df %>%
  filter(Fst > 0.4)

### LG15  NC_068435.1 18017912 - 19870544 tRNAs - bunch of genes
cum_df %>%
  filter(chr == "NC_068435.1" & midPos >=17500000 & midPos <=20500000)%>%
  ggplot() + 
  geom_point(aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  #geom_point(data = filter(lrrc9_df, Species == "Chum Even"), mapping = aes(x = cumPos, y = Fst), 
  #           color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  ylab("Chum") +
  geom_vline(xintercept=c(18017912,19870544),lty=2,color='#fe5a1d')

### LG23  NC_068443.1 10354488 - 10469110 kcnab1b
cum_df %>%
  filter(chr == "NC_068443.1" & midPos >=10000000 & midPos <=11000000)%>%
  ggplot() + 
  geom_point(aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  #geom_point(data = filter(lrrc9_df, Species == "Chum Even"), mapping = aes(x = cumPos, y = Fst), 
  #           color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  ylab("Chum") +
  geom_vline(xintercept=c(10354488,10469110),lty=2,color='#fe5a1d')

### LG28 NC_068448.1 850305 - 937353
cum_df %>%
  filter(chr == "NC_068448.1" & midPos >=750305 & midPos <=987353)%>%
  ggplot() + 
  geom_point(aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  #geom_point(data = filter(lrrc9_df, Species == "Chum Even"), mapping = aes(x = cumPos, y = Fst), 
  #           color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  ylab("Chum") +
  geom_vline(xintercept=c(850305,937353),lty=2,color='#fe5a1d')


### LG29  NC_068449.1 25440358 - 25545548 esrb
cum_df %>%
  filter(chr == "NC_068449.1" & midPos >=25000000 & midPos <=26545548)%>%
  ggplot() + 
  geom_point(aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  #geom_point(data = filter(lrrc9_df, Species == "Chum Even"), mapping = aes(x = cumPos, y = Fst), 
  #           color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  ylab("Chum") +
  geom_vline(xintercept=c(25440358,25545548),lty=2,color='#fe5a1d')

### LG35  NC_068455.1 28033205 - 28179022 lrrc9
cum_df %>%
  filter(chr == "NC_068455.1" & midPos >=27033205 & midPos <=28557254)%>%
  ggplot() + 
  geom_point(aes(x = midPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  #geom_point(data = filter(lrrc9_df, Species == "Chum Even"), mapping = aes(x = cumPos, y = Fst), 
  #           color = "brown1", alpha = 1, size = 1.5) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  ylab("Chum") +
  geom_vline(xintercept=c(28033205,28206418),lty=2,color='#fe5a1d')

esrb_df <- cum_df %>%
  filter(chr == "NC_068449.1" & midPos >=25000000 & midPos <=26545548)

lrrc9_df <- cum_df %>%
  filter(chr == "NC_068455.1" & midPos >=27033205 & midPos <=28557254)

er1_df <- cum_df %>%
  filter(chr == "NC_068428.1" & midPos >=23525841 & midPos <= 23552881)

Dummy_Leg <- cum_df %>% 
  filter((chr == "NC_068449.1" & midPos >=25000000 & midPos <=26545548) | (chr =="NC_068455.1" & midPos >=27033205 & midPos <=28557254) | (chr == "NC_068428.1" & midPos >=23525841 & midPos <= 23552881))%>%
  mutate(Color = ifelse(chr == "NC_068449.1" & midPos >=25000000 & midPos <=26545548,'esrb',
                        ifelse(chr == "NC_068455.1" & midPos >=27033205 & midPos <=28557254, 'lrrc9',
                               ifelse(chr == "NC_068428.1" & midPos >=23525841 & midPos <= 23552881,"esr1",NA))))%>%
  ggplot() + 
  geom_point(aes(x = cumPos, y = as.numeric(Fst), color = Color),
             alpha = 0.8, size = 1.3) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$LG[c(T,F)], 
                     breaks = axis_set$center[c(T,F)]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.05)) +
  ylab("Chum")+
  scale_color_manual(values = c('#4CBB17','dodgerblue2','#fe5a1d'),
                     labels = c('esr1','esrb','lrrc9'),
                     name="Peak")+
  theme(legend.position = 'right')

p_leg <- cowplot::get_legend(Dummy_Leg) 


##### Chum Manhattans ####
axis_set <- cum_df %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(center = mean(cumPos)) %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))

Chum_WG <- cum_df %>%
  mutate(chr_f = as.factor(chr),
         LG = as.numeric(chr_f))%>%
  ggplot() + 
  geom_point(aes(x = cumPos, y = as.numeric(Fst), color = as_factor(chr)),
             alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("gray50","black"),
                                  unique(length(axis_set$chr)))) +
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)), 
                     label = axis_set$LG[c(T,F)], 
                     breaks = axis_set$center[c(T,F)]) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0.2, 1, by = 0.2),
                     limits = c(0, 1.05)) +
  ylab("Chum") +
  theme(plot.margin = unit(c(0.1,0.1,0,0.1), "cm"),
        axis.text = element_text(angle = 0, size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  geom_point(data = lrrc9_df, mapping = aes(x = cumPos, y = as.numeric(Fst)), 
             color = "#fe5a1d", alpha = 1, size = 1.5)+
geom_point(data = esrb_df, mapping = aes(x = cumPos, y = as.numeric(Fst)), 
           color = "dodgerblue2", alpha = 1, size = 1.5)


DummyX <- ggplot() +
  labs(y = expression(italic(F['ST']))) + # Set the desired x-axis label
  theme(axis.line.x =  element_blank(),      # Keep the x-axis line
    axis.title.x =  element_blank(),     # Keep the x-axis title
    axis.text.x =  element_blank(),      # Keep the x-axis text (tick labels)
    #axis.title.y = element_blank(),    # Remove the y-axis title
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
pSpcAlgn <- cowplot::plot_grid(DummyX,
                   cowplot::plot_grid(pink_WG,
                                                          sockeye_WG,
                                                          Chum_WG,
                                                          coho_WG,
                                                          nrow=4,
                                                          labels = c("A", "B", "C", "D"),
                                                          label_fontfamily = "Times", 
                                                          label_colour = "black"),
                   p_leg,
                                       nrow=1,ncol=3,
                                       rel_widths = c(0.05,0.98,0.1))

png("./figures/FigureSX_SpeciesSpecificAlginment.png",height=8,width=11,res=300,units='in')
pSpcAlgn
dev.off()



