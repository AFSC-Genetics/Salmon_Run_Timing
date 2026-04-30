#### Summarize Protein Transcript Variants ####

SumGeneProtein <- function(
    Species,
    gene,
    Chrom ,
    #FstThreshold #considered omitting polymorphic sites that weren't different between groups
    Path2Fsts,
    Path2MAFearly,
    Path2MAFlate,
    Path2GeneSeq,
    Path2gff,
    OutputDir = "./results/mutations/",
    SaveTVrds = F,
    UseFstThreshold = T){
  
  #Preprocess Files
  # Create the CDS file for the species and gene of interest by subsetting the gff file
  #grep 'gene=lrrc9' /Volumes/PCCRC_v2/RunTiming_Salmon/Salmon_Run_Timing/data/raw/Genomes/PinkSalmonGenome/ncbi_dataset/data/GCF_021184085.1/genomic.gff > ./Pink_lrrc9.txt  
  gffcmd2pass <- paste("grep 'gene=",gene,";' ",Path2gff, " > ",OutputDir,gsub(pattern=" ",replacement="",Species),"_",gene,"_CDS.txt",sep="")  
  Path2CDS <- paste(OutputDir,gsub(pattern=" ",replacement="",Species),"_",gene,"_CDS.txt",sep="")
  system(gffcmd2pass,wait=T)  
  
  #Read in the gff file  
  gene_gff <- read.delim(Path2CDS,
                         header = F)%>%
    dplyr::rename('SeqID' = 'V1',
                  'Source' = 'V2',
                  'Type' = 'V3',
                  'Start' =  'V4',
                  'End' = 'V5',
                  'Score' = 'V6',
                  'Strand' = 'V7',
                  'Phase' =  'V8',
                  'Attrib' = 'V9')
  
  #Determine the start and end positions of the gene of interest
  GeneStartPos <- gene_gff %>% filter(Type=='gene') %>% select(Start) %>% unlist()
  GeneEndPos <- gene_gff %>% filter(Type=='gene') %>% select(End) %>% unlist()
  
  
  ### Step 1 Read in the Fst between early and late fish 
  
  Spcs_df <- read.delim2(Path2Fsts, header = T, sep = "\t", row.names = NULL) %>%
    mutate(Species = Species)%>%
    dplyr::select(region,chr,Nsites,Species)%>%
    dplyr::rename('chr' = 1,
                  'midPos' = 2,
                  'Fst' = 3)
  
  # find max of each chrom
  cumulate <- Spcs_df %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(max_bp = max(midPos)) %>%
    mutate(final_bp = max_bp + 5) %>%
    mutate(cum_bp = cumsum(final_bp))
  
  # create cumulative position column for all chromosomes
  cum_df <- NULL
  for(i in 1:length(unique(Spcs_df$chr))){
    chrom <- cumulate$chr[i]
    # for chr == 1 only, bc not adding to any chrom.
    if(i == 1){cum_temp <- Spcs_df %>%
      filter(chr == chrom) %>%
      mutate(cumPos = midPos)
    cum_df <- cum_temp
    }
    # for all other chromosomes, append
    if(i != 1){cumPos <- cumulate$cum_bp[i-1]
    cum_temp <- Spcs_df %>%
      filter(chr == chrom) %>%
      mutate(cumPos = midPos + cumPos)
    cum_df <- rbind(cum_df,cum_temp)
    }
  }
  
  
  #Subset to gene
  gene_df <- cum_df %>%
    filter(chr == Chrom & midPos >=GeneStartPos & midPos <=GeneEndPos) %>%
    mutate(Gene = gene)%>%
    mutate(Fst = as.numeric(Fst))%>%
    filter(Fst > 0)
  
  
  # I considered implementing a Fst threshold for summarizing, but decided to 
  # just list all polymorphic sites and number of sites where the minor allele differed
  if(UseFstThreshold == T){
    FstThreshold<-cum_df %>%
      filter(chr == Chrom)%>%
      mutate(Fst = as.numeric(Fst))%>%
      mutate(Fst = case_when(
        Fst < 0 ~ 0,
        .default = Fst))%>%
      summarize(meanFst = mean(as.numeric(Fst),na.rm=T), 
                sdFst=sd(as.numeric(Fst),na.rm=T),
                FstThresh = meanFst+2*sdFst)
    
    gene_df <- gene_df %>%
      filter(Fst >= FstThreshold$FstThresh)
    
  }
  
  
  
  #### Step 2 read in the MAF for the early and late fish
  MAFearly <- read_delim(Path2MAFearly)
  MAFlate <- read_delim(Path2MAFlate)
  
  MAFpheno <- MAFearly %>%
    full_join(.,MAFlate,by=c('chromo','position','ref','anc'))%>%
    arrange(position)
  
  #Bind the polymorphic sites from the study to the MAF data 
  gene_df <- gene_df %>%
    dplyr::rename('chromo' = 'chr',
                  'position' = 'midPos') %>%
    left_join(.,MAFpheno,by=c('chromo','position')) #this was a full join
  
  #Make sure we are just looking at the gene of interest 
  gene_polymorphic <- gene_df %>%
    filter(Gene == gene) 
  
  #Read in the reference sequence of the gene downloaded from NCBI (https://www.ncbi.nlm.nih.gov/nuccore/NC_060182.1?report=genbank&from=72710811&to=72735126)
  GeneSeq <- readLines(Path2GeneSeq %>% unlist() %>% file.path())%>%
    lapply(.,str_split,"") %>%
    unlist()
  
  
  ### Step 3 get the CDS info from the annotation file 
  # grep 'gene=lrrc9' /Volumes/PCCRC_v2/RunTiming_Salmon/Salmon_Run_Timing/data/raw/Genomes/PinkSalmonGenome/ncbi_dataset/data/GCF_021184085.1/genomic.gff > ./Pink_lrrc9.txt
  
  #Check to see if there is a CDS
  Types <- gene_gff %>%
    select(Type) %>%
    filter(!Type%in%c('gene','exon'))%>%
    unlist() %>%
    as.vector()
  
  if (any(grepl('CDS|cds',Types))){
    #Get the CDS info for transcribing the protein
    gene_CDS <-  gene_gff %>%  
      filter(Type=="CDS") %>%
      mutate(Transcript_Variant =  str_extract(Attrib, "(?<=protein_id=)XP_[^;]+"))
    
    Isoform_Bases <- lapply(1:length(unique(gene_CDS$Transcript_Variant)),function(x)
      gene_CDS %>%
        filter(Transcript_Variant == unique(gene_CDS$Transcript_Variant)[x])%>%
        rowwise()%>%
        mutate(Bases = paste(seq(from=Start,to=End,by=1),collapse=","))%>%
        select(Bases)%>%
        unlist()%>% paste(.,collapse=",")%>%
        lapply(.,str_split,",")%>%
        unlist())
    
    #### Step 4 Bind the Gene Sequence to the polymorphic sites
    
    #Bind the full gene sequence with allele freq data from our collections
    GeneSeq_df <- data.frame(position=seq(from = gene_gff %>% filter(Type=='gene') %>% select(Start) %>% unlist(),
                                          to = gene_gff %>% filter(Type=='gene') %>% select(End) %>% unlist(),
                                          by=1),
                             Base = GeneSeq[-1])%>%
      left_join(.,gene_polymorphic,by=c('position'))
    
    #Now lets, for each transcript variant, identify which bases are in each CDS isoform
    for (TV in unique(gene_CDS$Transcript_Variant)){
      TVc <- which(unique(gene_CDS$Transcript_Variant) %in% TV) # counter
      GeneSeq_df<- GeneSeq_df %>%
        mutate(!!TV := ifelse(position %in% as.numeric(unlist(Isoform_Bases[[TVc]])),1,NA))
    }
    
    # How many sites in gene, cds, polymorphic in cds?
    TVsumSites <- GeneSeq_df %>%
      select(position,Base,major.x,minor.x,unique(gene_CDS$Transcript_Variant)) %>%
      reshape2::melt(id.vars=c('position','Base','major.x','minor.x'))%>%
      mutate(Polymorphic = ifelse(is.na(major.x),NA,1))%>%
      group_by(variable) %>%
      summarize(nSitesGene = n(),
                nSitesCDS = sum(value,na.rm=T),
                nSitesGenePoly = sum(Polymorphic,na.rm=T))
    
    #of the polymorphic sites in the gene, how many are polymorphic in the CDS
    CDSsumSites <- GeneSeq_df %>%
      select(position,Base,major.x,minor.x,unique(gene_CDS$Transcript_Variant)) %>%
      reshape2::melt(id.vars=c('position','Base','major.x','minor.x'))%>%
      mutate(Polymorphic = ifelse(is.na(major.x),NA,1))%>%
      filter(value == 1) %>% #only select those cites in the CDS
      group_by(variable) %>%
      summarize(nSitesCDS = sum(value,na.rm=T),
                nSitesCDSPoly = sum(Polymorphic,na.rm=T))
    
    GeneSumSites <- TVsumSites %>%
      left_join(.,CDSsumSites,by=c('variable','nSitesCDS'))
    
    if(all(GeneSumSites$nSitesCDSPoly == 0)){
      NoPolyCDSres <- cbind(data.frame(Species=Species,
                                       Gene=gene),
                            GeneSumSites)%>%
        as_tibble()%>%
        dplyr::rename("TV"='variable')%>%
        mutate(Type='CDS')%>%
        return(NoPolyCDSres)
    } else {
      
      #Plug in some degenerate bases for the Consensus base
      gene_CDS_df <- GeneSeq_df %>%
        select(position,Base,major.x,minor.x,unique(gene_CDS$Transcript_Variant)) %>%
        reshape2::melt(id.vars=c('position','Base','major.x','minor.x')) %>%
        filter(value == 1) %>% #only do it for the CDS
        mutate(CnsnsBase = ifelse(is.na(major.x),Base,paste(major.x,minor.x,sep="/")))%>%
        mutate(DgnrtBase = case_match(CnsnsBase,
                                      c("A/G", "G/A") ~ "R",
                                      c("C/T", "T/C") ~ "Y",
                                      c("A/C", "C/A") ~ "M",
                                      c("G/T", "T/G") ~ "K",
                                      c("G/C", "C/G") ~ "S",
                                      c("A/T", "T/A") ~ "W",
                                      .default = CnsnsBase # Keep original if no match
        ))
      
      #What if we wanted to know the sequence for the Early and late allele
      
      PolyGenBs <- gene_CDS_df %>%
        filter(grepl(pattern="/",CnsnsBase))%>%
        left_join(.,gene_df %>% select(position,Fst),by="position")
      
      #what if there are no mutations in the coding region?
      
      #Early
      # AllelesVT <-list()
      # cnt<-1
      # for (TV in unique(gene_CDS$Transcript_Variant)){
      # EAlist<-list()
      # LAlist<-list()
      # 
      # for (r in 1:nrow(PolyGenBs)){
      #   PosTemp <- PolyGenBs %>%
      #     filter(variable == TV) %>%
      #     select(position)%>% unlist() %>%
      #     {.[r]} 
      #   
      #   EAlist[[r]] <- MAFearly %>%
      #     filter(position==PosTemp) %>%
      #     mutate(EA = ifelse(knownEM<0.5,major,minor))%>%
      #     mutate(variable = TV) %>%
      #     select(chromo,position,variable,EA)
      # 
      #   LAlist[[r]] <- MAFlate %>%
      #     filter(position==PosTemp) %>%
      #     mutate(LA = ifelse(knownEM<0.5,major,minor))%>%
      #     mutate(variable = TV) %>%
      #     select(chromo,position,variable,LA)
      # }
      # 
      # AllelesVT[[cnt]] <- full_join(do.call(rbind,EAlist) %>% select(chromo,position,EA,variable),
      #                                                       do.call(rbind,LAlist)%>% select(chromo,position,LA,variable),
      #                                                       by=c('chromo','position','variable'))
      #   cnt<-1+cnt
      # }
      
      # 1. Prepare Early and Late alleles globally (vectorized)
      EA_df <- MAFearly %>%
        mutate(EA = if_else(knownEM < 0.5, major, minor)) %>%
        select(chromo, position, EA)
      
      LA_df <- MAFlate %>%
        mutate(LA = if_else(knownEM < 0.5, major, minor)) %>%
        select(chromo, position, LA)
      
      # 2. Join everything at once
      # This replaces the entire nested loop structure
      AllelesVT_Final <- PolyGenBs %>%
        mutate(chromo = Chrom)%>%
        #dplyr::rename('TV' = 'variable') %>% # ensure column name matches your logic
        left_join(EA_df, by = c("chromo", "position")) %>%
        left_join(LA_df, by = c("chromo", "position"))
      
      # 3. If you specifically need it as a list split by Transcript_Variant:
      AllelesVT <- split(AllelesVT_Final, AllelesVT_Final$variable)
      
      
      gene_CDS_df2 <- gene_CDS_df %>%
        left_join(.,do.call(rbind,AllelesVT) %>% select(position,EA,LA,variable),by=c('position','variable'))%>%
        left_join(.,gene_polymorphic %>% select(position,Fst),by=c('position')) %>%
        mutate(EA=ifelse(is.na(EA),CnsnsBase,EA),
               LA=ifelse(is.na(LA),CnsnsBase,LA))
      
      MutsinCDS <-gene_CDS_df2 %>%
        mutate(Diffs = ifelse(EA==LA,0,1))%>%
        group_by(variable) %>%
        summarize(nDiffsELminor=sum(Diffs))
      
      
      #Now lets write out the fasta files for each transcipt
      TVseq<-list()
      cnt<-1
      for (al in c('DgnrtBase','EA','LA')){
        for (TV in unique(gene_CDS$Transcript_Variant)){
          TVseq[[cnt]] <- gene_CDS_df2 %>%
            filter(variable == TV) %>%
            arrange(position) %>%
            select(all_of(al)) %>%
            unlist() %>%as.vector() %>% paste(.,collapse="")
          
          names(TVseq[[cnt]])<-paste(al,TV,sep="-")
          
          cnt<-cnt+1
        }
      }
      
      # names(TVseq) <- crossing(c('DgnrtBase','EA','LA'),unique(gene_CDS$Transcript_Variant))%>%
      #   dplyr::rename('Allele'=1,
      #          'TV' = 2)%>%
      #   rowwise() %>%
      #   mutate(name = paste(Allele,TV,sep="-"))%>%
      #   select(name)%>%
      #   unlist()
      
      #which strand is the protein coded from?
      ForRstrand <- ifelse(gene_CDS$Strand %>% unique()=="+","F",
                           ifelse(gene_CDS$Strand %>% unique()=="-","R",
                                  stop(paste("Check the strand for gene ", gene," and species", Species,sep=""))))
      
      #Lets now convert to a protein product
      TVpro<-list()
      for (x in 1:length(TVseq)){
        TVpro[[x]] <-  TVseq[[x]] %>%
          unlist() %>%
          str_split(.,"") %>%
          unlist()%>%
          seqinr::translate(.,frame=0,sens=ForRstrand,numcode=1,ambiguous=T)%>%
          paste(.,collapse="")
        
        names(TVpro[[x]])<-names(TVseq[[x]])
        
      }
      
      # names(TVpro) <- crossing(c('DgnrtBase','EA','LA'),unique(gene_CDS$Transcript_Variant))%>%
      #   dplyr::rename('Allele'=1,
      #          'TV' = 2)%>%
      #   rowwise() %>%
      #   mutate(name = paste(Allele,TV,sep="-"))%>%
      #   select(name)%>%
      #   unlist()
      
      
      TVlst <- list()
      TVel <- lapply(TVpro,names)%>%unlist() %>% grep(pattern="EA|LA",value = T)
      
      TVpro_bs <- TVpro[which(lapply(TVpro,names)%>%unlist()%in%TVel)]%>%
        str_split(.,"")
      names(TVpro_bs) <- TVel
      
      
      #Save them if requested
      if(SaveTVrds == T){
        if(!dir.exists(paste(OutputDir,"TranscriptVariants",sep=""))){
          dir.create(file.path(paste(OutputDir,"TranscriptVariants",sep="")))
        }
        
        saveRDS(file=paste(OutputDir,"TranscriptVariants/",Species,
                           "_",gene,"_",Chrom,".RDS",sep=""),
                object=TVpro[which(lapply(TVpro,names)%>%unlist()%in%TVel)])
      }
      
      cnt<-1 
      for (TV in unique(gene_CDS$Transcript_Variant)){
        TVlst[[cnt]] <- TVpro_bs[grep(pattern=TV,x=names(TVpro_bs))] %>%
          {lapply((1:2),function(x) data.frame(Allele=unlist(lapply(str_split(names(.[x]),"-"),"[[",1)),
                                               TV = unlist(lapply(str_split(names(.[x]),"-"),"[[",2)),
                                               Seq = unlist(.[x]),
                                               pos = seq(1,length(unlist(.[x])),by=1)))}%>%
          {left_join(.[[1]],.[[2]],by=c('pos','TV'))}%>%
          dplyr::rename('EarlyP' = 'Seq.x' ,
                        'LateP' = 'Seq.y' )%>%
          select(pos,TV,EarlyP,LateP)%>%
          mutate(Diffs = ifelse(EarlyP==LateP,0,1))
        
        cnt<-1+cnt
        
      }
      
      ProSum <- TVlst %>%
        do.call(rbind,.)%>%
        group_by(TV) %>%
        # lapply(1:length(TVlst),function(x)
        #   TVlst[[x]] %>%
        #     group_by(TV) %>%
        summarize(nAA = n(),
                  nNYmuts = sum(Diffs),
                  posDiffs = paste(pos[Diffs == 1], collapse = "-"),
                  muts = paste(paste0(EarlyP[Diffs == 1], "-", LateP[Diffs == 1]), collapse = ","))
      
      #FstDiffs?
      gene_CDS_df2_with_codons <- gene_CDS_df2 %>%
        group_by(variable) %>%
        mutate(codon = {
          # Calculate codons for this specific group
          n_positions <- n()
          n_codons <- n_positions / 3
          
          if (ForRstrand == "F") {
            rep(seq_len(n_codons), each = 3)
          } else if (ForRstrand == "R") {
            rep(seq(from = n_codons, to = 1, by = -1), each = 3)
          } else {
            NA
          }
        }) %>%
        ungroup()
      
      CodonMuts <- gene_CDS_df2_with_codons %>%
        filter(grepl(pattern="/",CnsnsBase))
      
      # Summary of gene
      SumTib <- ProSum %>% #number of amino acids, number non-synonymous mutations, positions of the mutations, mutations
        left_join(GeneSumSites %>% dplyr::rename("TV" = variable),by="TV")%>%
        left_join(., MutsinCDS%>% dplyr::rename("TV" = variable),by="TV")%>%
        mutate(Species=Species,
               Gene=gene)%>%
        select(Species,Gene,TV,nSitesGene,nSitesGenePoly,nSitesCDS,nSitesCDSPoly,nDiffsELminor,nAA,nNYmuts,posDiffs,muts)%>%
        dplyr::rename('nNSynmuts' = 'nNYmuts' )
      
      SumTibSNP <-list()
      for (st in 1:nrow(SumTib)){
        SumTibSNP[[st]] <- CodonMuts %>%
          filter(variable == SumTib[st,]$TV)%>%
          filter(codon %in% unlist(str_split(SumTib[st,]$posDiffs,"-")))%>%
          dplyr::rename('posDiffs' = 'codon',
                        'TV'='variable')%>%
          select(TV,posDiffs,position,Fst)%>%
          arrange(posDiffs)%>% #this arranges the positions in the codon that differ so that we can bind correctly later. 
          group_by(TV) %>% 
          summarize(posDiffs = paste(posDiffs,collapse="-"),
                    positions = paste(position,collapse="-"),
                    Fst = paste(Fst,collapse="-"))
      }
      SumTibSNP <- do.call(rbind,SumTibSNP)
      
      SumTib <- SumTib %>%
        left_join(.,SumTibSNP,by=c('TV','posDiffs'))%>%
        mutate(Type = "CDS")%>%
        relocate(Type,.after='Gene')
      
      return(SumTib)
    } # over polymorphic sites in CDS
  }else{ #if there is a CDS in the gff
    SumTib <- data.frame(Species=Species,
                         Gene = gene,
                         Type = paste(Types,collapse=""))
    return(SumTib)
  }
}

