# LOCAL SCORE
# input is a table with p-vals from a FET based on counts data between groups 1 and 2

#title: "R script for computing local scores"
# POOLPARTY PAPER LINK: https://onlinelibrary.wiley.com/doi/full/10.1111/mec.16810

### **Description**: 
## Detecting genomic outliers in a sequence of correlated scores, from Fariello et al 2017.
#  LOCAL SCORE PAPER LINK: https://onlinelibrary.wiley.com/doi/10.1111/mec.14141 

### Import libraries ##########################################################
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(tidyverse)

### Inputs ####################################################################

# choose xi (may have to be adjusted)
filename="pink-chum_minInd0.3_MM3"

# fishers exact test output
fet_df <- fread(paste0("results/localscore/",filename,"_postFET.tsv"), head=T)

### Functions - Local Score ############################################# ##################################
# computation of the autocorrelation
autocor=function(x){abs(cor(x[-1],x[-length(x)]))} 

### computation of the lindley process from scores
lindley=function(scores){
  L=length(scores)
  sl=rep(0,L+1)
  for (i in 1:L){
    sl[i+1]= max(0, (sl[i]+scores[i]))
  }
  return(sl[-1])
}

# computation of the local score from a lindley process and the region under 
# the maximum peak
scorelocal=function(lind){
  nC=max(lind[,1])
  sloc=matrix(0, nrow=nC, ncol=3)
  colnames(sloc)=c('chr','beg','end')
  for (i in 1:nC){
    tmp=which(lind[,1]==i)		
    list=lind[tmp,3]
    M_loc=which.max(list)
    if(length(which(list[1:M_loc]==0))>0){
      m_loc=max(which(list[1:M_loc]==0))+1			
    }else{m_loc=1}
    sloc[i,]=c(i,lind[tmp[m_loc],2],lind[tmp[M_loc],2])	
  }
  return(data.frame(zone=seq(1,nrow(sloc)),sloc))
}

#computation of the significance threshold if the distribution of the p-values is uniform
thresUnif=function(L, cor, xi, alpha = 0.05){
  coefs=list(list('a'=c(-5.5,6.76,-5.66,-2.51),'b'=c(-1.22,3.17,-1.99)),
             list('a'=c(2.47,-4.16,-1.82,-4.58),'b'=c(0.37,2.14,-2.35)),
             list('a'=c(2.04,-5.76,1.04,-6.95),'b'=c(2.55,-0.02,-2.31)),
             list('a'=c(0.22,-4.08,1.16,-9.16),'b'=c(3.45,-0.98,-2.33))
  )
  cors=c(cor^3,cor^2,cor,1)
  if (missing(xi) | !(xi %in% 1:4)){
    print('xi should be 1, 2, 3 or 4') 
  }else{
    a=log(L)+ coefs[[xi]]$a %*%cors
    b=coefs[[xi]]$b %*%cors[-1]
    #then compute the threshold:
    thres = ( log(-log(1-alpha)) - a ) / b
    return(thres)
  }
}

# computation of the significative regions from a lindley process given a significance threshold
sig_sl=function(lind,pos, th){
  zones=c(0,0,0)	
  list=lind
  auxpos=pos
  while(max(list)>=th){
    M_loc=which.max(list)
    if(length(which(list[1:M_loc]==0))==0){ #the peak is at the beginning of the chrom
      m_loc=1
      zones=rbind(zones, c(auxpos[m_loc],auxpos[M_loc],max(list)))
      tmp=which.min[which(list[M_loc+1:length(list)]==0)] #first 0 score after peak
      list=list[tmp:length(list)]
      auxpos=pos[tmp:length(list)]
    }else{	
      m_loc=max(which(list[1:M_loc]==0))			
      max=max(list)
      zones=rbind(zones, c(auxpos[m_loc+1],auxpos[M_loc],max))
      tmp=which(list[M_loc:length(list)]==0) #first 0 score after peak
      if (length(tmp)>0){
        auxpos=auxpos[c(1:m_loc,(min(tmp)+M_loc):length(list))]
        list=list[c(1:m_loc, (min(tmp)+M_loc):length(list))]
      }else{ #if the peak is at the end of the chromosome
        auxpos=auxpos[1:m_loc]
        list=list[1:m_loc]
      }				
    }
  }
  zones=matrix(zones, ncol=3)
  zones=data.table(beg=zones[,1],end=zones[,2],peak=zones[,3])
  if (nrow(zones)>1){zones=zones[-1,]}
  return(zones)
}

### Estimation of Gumble coefficients.

# Estimation of the coefficients of the Gumbel distributions
lineGumb=function(x){
  x1tmp=x
  if (length(table(x1tmp)) >  5){
    Fn = ecdf(x1tmp)
    Fnx=Fn(x1tmp)
    filtre = ( Fnx < 1 - 10**(-1) )  & ( x1tmp > 1 ) #sinon pbs numq
    lm0 = lm(log(-log(Fnx[filtre])) ~ x1tmp[filtre])$coefficients
  }else{lm0=c(0,0)}  
  return(lm0)
}

# Estimation of the coefficients of the polynomes to compute the Gumbel coefficients 
# depending on the length of the chromosomes and the chromosome autocorrelation
coefsGumb=function(mydata, Ls=seq(10000,80000,10000), nSeq=5000, degMod=2){
  bins=seq(0.025,0.975,0.05) 
  coefs=array(0, dim=c(length(bins),5, length(Ls)))  
  coefs[,5,]=matrix(rep(Ls,length(bins)), ncol=length(Ls),byrow=TRUE)
  as=matrix(0, ncol=length(Ls), nrow=length(bins))
  bs=matrix(0, ncol=length(Ls), nrow=length(bins))
  xx=seq(0,max(mydata$lindley), 1)
  for (j in 1:length(Ls)){
    coefs[,1,j]=bins  
    len=Ls[j]
    tmp=sample(nrow(mydata)-len,nSeq,replace=F) #samples the sequences'beginnings
    DTL=data.table(LocScore = vector('numeric', length=nSeq),cors= vector('numeric', length=nSeq)) 
    for (l in 1:nSeq){
      DTL[l]=mydata[seq(tmp[l],length.out=len),.(max(lindley(score)),autocor(pval))]
    }
    DTL[,bin:=which.min(abs(bins-cors)),cors]
    binNE=sort((unique(DTL$bin)))
    coefs[binNE,4,j]=table(DTL$bin)
    DTCoef=DTL[,.(lineGumb(LocScore)),bin]
    coefs[unique(DTL$bin),c(2,3),j]=matrix(DTCoef$V1, ncol=2, byrow=T)
    ys=coefs[,c(2,3),j]%*%rbind(rep(1,length(xx)), xx)
    pdf(paste('GumbelinesAllxi',xi,'M',len,'.pdf', sep=''))
    for (i in binNE){
      Fn=ecdf(DTL[bin==i,LocScore]) 
      tmp=DTL[bin==i,.(LocScore,Fn=Fn(LocScore)),][(Fn< 1 - 10**(-5)) & (LocScore > 1)]
      if(i==min(binNE)){
        plot(tmp$LocScore, log(-log(tmp$Fn)), col=i,main=paste("M = ",len,sep=""), xlim=quantile(DTL$LocScore,c(0,0.9)),ylim=range(ys[binNE,xx<quantile(DTL$LocScore,0.8)], na.rm=T), pch=20, xlab='Local Score', ylab='log(-log(Fn(LS)))')
      }else{points(tmp$LocScore, log(-log(tmp$Fn)), col=i, pch=20)
      }
      if (sum(ys[i,]!=0)){lines(xx, ys[i,], col=i, lwd=1.5)}
    }
    legend("topright",legend=bins[binNE], pch=16,col=seq(1:length(binNE)), title='rho')
    dev.off()
    as[,j]=as.numeric(coefs[,2,j])
    bs[,j]=as.numeric(coefs[,3,j])		
  }
  
  auxWhich=which(coefs[,4,]> 150)
  rhos=coefs[,1,][auxWhich]
  auxAs=as[auxWhich]-log(coefs[,5,][auxWhich])
  auxBs=bs[auxWhich]
  
  fitA=lm((auxAs)~I(rhos))  
  fitB=lm((auxBs)~I(rhos))  
  
  pdf('FitAandB.pdf')
  par(mfrow=c(1,2))
  plot(rhos, auxAs, pch=16, col=coefs[,5,][auxWhich]/10000, xlab='rho', ylab='a-log(M)')
  xslm=seq(0,1,0.01)
  lines(xslm, fitA$coefficients%*%rbind(rep(1,length(xslm)), xslm))
  plot(rhos, auxBs, pch=16, col=coefs[,5,][auxWhich]/10000, xlab='rho', ylab='b')
  lines(xslm, fitB$coefficients%*%rbind(rep(1,length(xslm)), xslm))
  legend(min(rhos), max(auxBs), legend=sort(unique(coefs[,5,][auxWhich])), col=unique(coefs[,5,][auxWhich]/10000), pch=16, title='M values')
  
  dev.off()
  
  return(list(aCoef=fitA$coefficients, bCoef=fitB$coefficients))
}

# Computation of the significance threshold given the computed polynomes for computing
# the Gumbel coefficients depending of length and autocorrelation
threshold=function(L, cor, aCoef, bCoef, alpha = 0.05){
  degA=length(aCoef)
  degB=length(bCoef)
  a=log(L)
  b=0
  for (i in 1:degA){
    a=a+aCoef[i]*cor^(i-1)  
  }
  for (i in 1:degB){
    b=b+bCoef[i]*cor^(i-1)  
  }
  #then compute the threshold:
  thres = ( log(-log(1-alpha)) - a ) / b
  return(thres)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# END OF LOCALSCOREFUNCTIONS.R



### POOLPARTY2 Code ###########################################################

# Read in data
#chromosome number | position | score
#The first line of the table should include the names: chr pos pval
setkey(fet_df, chr)
fet_df$pval[fet_df$pval==0]=1e-16
Nchr=length(fet_df[,unique(chr)])
head(fet_df)

### Computation of absolute position in the genome (posT)
#This is useful for doing genomewide plots. 
chrInfo=fet_df[,.(L=.N,cor=autocor(pval)),chr]
setkey(chrInfo,chr)
data.table(chr=fet_df[,unique(chr),], S=cumsum(c(0,chrInfo$L[-Nchr]))) %>% 
  setkey(.,chr) %>% fet_df[.,posT:=pos+S]

# See p-value distributions

pname= paste0("figures/localscore/parameterization/",filename,"_pvaldist.pdf")
pdf(pname)
qplot(pval, data=fet_df, geom='histogram', binwidth=0.1, main='P-values histogram')
invisible(dev.off())

bname= paste0("figures/localscore/parameterization/",filename,"_blog10.pdf")
pdf(bname)
qplot(-log10(pval), data=fet_df, geom='histogram', binwidth=0.1, main='-log10(p-values) histogram')
invisible(dev.off())

meanP= mean(fet_df$pval)
stdP= sd(fet_df$pval)

#Kolmogorv Smirnov Test for uniformity

minP=min(fet_df$pval)
maxP=max(fet_df$pval)

KST=ks.test(unique(fet_df$pval), "punif",minP,maxP)[[2]]

meanl10<- mean(-log10(fet_df$pval))
print(paste0("Mean of -log10p is ", round(meanl10,3)))

## Computation of the score and the Lindley Process


### Xi vs Quantile Relationship ###############################################
print("The Xi (or tuning) parameter controls the degree of smoothing the underlying distribution of significance: ")
print("It represents a tradeoff between precision and power: ")
print("Lower value (less smoothing) provides more power but less precision (larger islands of divergence; more false positives);")
print("Larger values (more smoothing) provide greater precision, but less power (more false negatives).")
print("The quantiles of -log10p (potential choices of Xi) are (also printed to file): ")
print(c(quantile(-log10(fet_df$pval), probs = seq(0.6, 0.95, 0.05)),quantile(-log10(fet_df$pval), probs = seq(0.96, 0.99, 0.01))))

print("Which produce differences of (only negative values will work!): ")
print(meanl10-c(quantile(-log10(fet_df$pval), probs = seq(0.6, 0.95, 0.05)),quantile(-log10(fet_df$pval), probs = seq(0.96, 0.99, 0.01))))

parameterization <- paste0("results/localscore/parameterization/",filename,"_tuning_parameter_xi_percentiles.txt")
write.table("The quantiles of -log10p (potential choices of Xi) are: ",
            file=parameterization, col.names=F, row.names=F, quote=F, append=F,sep="\t")
write.table(rbind(as.integer((c(seq(0.6, 0.95, 0.05),seq(0.96, 0.99, 0.01))*100)),
                  c(quantile(-log10(fet_df$pval), probs = seq(0.6, 0.95, 0.05)),quantile(-log10(fet_df$pval), probs = seq(0.96, 0.99, 0.01)))),
            file=parameterization, col.names=F,row.names=F,quote=F,append=T,sep="\t")
write.table("Which produce differences of (only negative values will work!): ",
            file=parameterization, col.names=F,row.names=F,quote=F,append=T,sep="\t")
write.table(rbind(as.integer((c(seq(0.6, 0.95, 0.05),seq(0.96, 0.99, 0.01))*100)),
                  meanl10-c(quantile(-log10(fet_df$pval), probs = seq(0.6, 0.95, 0.05)),quantile(-log10(fet_df$pval), probs = seq(0.96, 0.99, 0.01)))),
            file=parameterization, col.names=F,row.names=F,quote=F,append=T,sep="\t")
write.table("The smallest possible quantile of Xi is: ",
            file=parameterization,col.names=F,row.names=F,quote=F,append=T,sep="\t")
write.table(rbind((min(which(meanl10-quantile(-log10(fet_df$pval), probs = seq(0.001, 0.999, 0.001))<0))/10),
                  round(quantile(-log10(fet_df$pval), probs = seq(0.001, 0.999, 0.001)),4)[min(which(meanl10-quantile(-log10(fet_df$pval), probs = seq(0.001, 0.999, 0.001))<0))]),
            file=parameterization, col.names=F,row.names=F,quote=F,append=T,sep="\t")

# CHECK xi choice
# To choose the appropriate threshold ($\xi = 1, \dots, 4$) we look at the distribution of −log10(p−value). 
#   Then the score function will be $−log10(p−value) − \xi$. 

# We should verify that the mean of the score is negative. 
#It should be according to the chosen value of xi.
xi=3
fet_df[,score:= -log10(pval)-xi]
fet_df[,lindley:=lindley(score),chr]

# The score mean must be negative to proceed!
if(mean(fet_df$score) >= 0){
  print(paste0("Error: xi needs to be reset so the mean < 0! Mean = ",mean(fet_df$score)))
}else{
  print(paste0("Passed: Mean of xi is negative. Mean = ",mean(fet_df$score)))
}


### Non-Uniform #############
# re-sampling strategy using Gumbel's coefficients
# Gumbel distribution can be applied to represent the distribution of maxima

# If you see this error:
#  "Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) :0 (non-NA) cases"
# It usually indicates a sampling error; try re-running the script.

# Calculating non-uniform coefficients
minSNPsPerChrom=min(table(fet_df$chr))
RoundedMinSNPsPerChrom=max(seq(10000,80000,10000)[seq(10000,80000,10000)<=minSNPsPerChrom])
RoundedMinSNPsPerChrom

coefsG=coefsGumb(fet_df, Ls=seq(10000,RoundedMinSNPsPerChrom,5000), nSeq=5000)

# Determining significance thresholds and writing files
chrInfo[,thG05:=threshold(L, cor, coefsG$aCoef, coefsG$bCoef,0.05),]
chrInfo[,thG01:=threshold(L, cor, coefsG$aCoef, coefsG$bCoef,0.01),]
chrInfo[,thG001:=threshold(L, cor, coefsG$aCoef, coefsG$bCoef,0.001),]
chrInfo[,thG075:=threshold(L, cor, coefsG$aCoef, coefsG$bCoef,0.075),]

chrInfo
# save chromosome level significance thresholds
write.table(chrInfo, file = paste0("results/localscore/parameterization/",filename,"_xi",xi,"_chrs_signif.txt"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

fet_df=fet_df[chrInfo]
fet_df=fet_df[thG05>0&thG01>0&thG001>0&thG075>0]

sigZones05=fet_df[,sig_sl(lindley, pos, unique(thG05)),chr]
sigZones01=fet_df[,sig_sl(lindley, pos, unique(thG01)),chr]
sigZones075=fet_df[,sig_sl(lindley, pos, unique(thG075)),chr]
sigZones001=fet_df[,sig_sl(lindley, pos, unique(thG001)),chr]

sig05name = paste0("results/localscore/parameterization/",filename,"_xi",xi, "_05sig.txt")
sig01name = paste0("results/localscore/parameterization/",filename,"_xi",xi, "_01sig.txt")
sig001name = paste0("results/localscore/parameterization/",filename,"_xi",xi, "_001sig.txt")
sig075name = paste0("results/localscore/parameterization/",filename,"_xi",xi, "_075sig.txt")

ind=which(sigZones075[,peak]>0)
if (nrow(sigZones075) >0 ) {
  write.table(sigZones075[ind,],file=sig075name,col.names=T,row.names=F,quote=F,sep="\t")
}
ind=which(sigZones05[,peak]>0)
if (nrow(sigZones05) >0 ) {
  write.table(sigZones05[ind,],file=sig05name,col.names=T,row.names=F,quote=F,sep="\t")
}
ind=which(sigZones01[,peak]>0)
if (nrow(sigZones01) >0 ) {
  write.table(sigZones01[ind,],file=sig01name,col.names=T,row.names=F,quote=F,sep="\t")
}
ind=which(sigZones001[,peak]>0)
if (nrow(sigZones001) >0 ) {
  write.table(sigZones001[ind,],file=sig001name,col.names=T,row.names=F,quote=F,sep="\t")
}

#Average of significance thresholds --> write to table
cal05 <- mean(fet_df$thG05)
cal01 <- mean(fet_df$thG01)
cal075 <- mean(fet_df$thG075)
cal001 <- mean(fet_df$thG001)

avname = paste0("results/localscore/parameterization/",filename,"_xi",xi, "_mean_sig.txt")
AVsig=rbind(c(0.001,0.01,0.05,0.075), c(cal001,cal01,cal05,cal075))
write.table(AVsig,file=avname,col.names=F,row.names=F,quote=F,sep="\t")

#Chromosome specific thresholds 
fet_dfchr <- fet_df[,unique(fet_df$chr)]
fet_dfchr<- fet_df[!duplicated(fet_df$chr), ]

cal05chr <- fet_dfchr$thG05
cal01chr <- fet_dfchr$thG01
cal075chr <- fet_dfchr$thG075
cal001chr <- fet_dfchr$thG001

chrname = paste0("results/localscore/parameterization/",filename,"_xi",xi,"_chr_sig.txt")
CHsig<-as.data.frame(cbind(fet_dfchr$chr,cal001chr,cal01chr,cal05chr,cal075chr))
colnames(CHsig) <- c("Chr", "sig.001", "sig.01", "sig.05", "sig.075")
write.table(CHsig,file=chrname,col.names=T,row.names=F,quote=F,sep="\t")

lsname = paste0("results/localscore/parameterization/",filename,"_xi",xi, ".ls")
lsOUT<- as.data.frame(cbind(fet_df$chr, fet_df$pos, fet_df$snp, fet_df$lindley))
write.table(lsOUT,file=lsname,col.names=F,row.names=F,quote=F,sep="\t")

# Plotting chr-by-chr
pdf(paste0("figures/localscore/",filename,'_ScoreLocalAll_xi',xi,'_nonuniform.pdf'))
par(mfrow=c(4,1), mar=c(5,2,1,1))
for (g in chrInfo$chr){
  plot(fet_df[chr==g,pos],fet_df[chr==g,lindley],
       type='l',
       xlab=paste("position chr",g,sep=' '),
       ylab="Score Local",
       #xlim=c(0,max()),
       ylim=c(0,max(fet_df[chr==g,max(lindley)],fet_df[chr==g, unique(thG01)])))
  abline(h=fet_df[chr==g, unique(thG05)], col='dodgerblue3')
  abline(h=fet_df[chr==g, unique(thG01)], lty=3, col='firebrick4')
  abline(v=sigZones01[chr==g, beg], col='pink3', lty=3)
  abline(v=sigZones01[chr==g, end], col='pink4', lty=2)
}
dev.off()


### Uniform Distribution ######################################################
# # Compute significance threshold for each chromosome
# 
# ## Uniform distribution of p-values
# 
# #If the distribution of the p-values is uniform and if $\xi=$ is 1, 2, 3 or 4, 
# #it is possible to compute the significancy thresholds for each chromosome directely, 
# #given the length and the autocorrelation of the chromosome .
# if(KST > 0.5){
#   chrInfo[,th:=thresUnif(L, cor, 1),chr]
#   fet_df=fet_df[chrInfo]
#   
#   sigZones=fet_df[chrInfo,sig_sl(lindley, pos, unique(th)),chr]
#   
#   pdf(paste0("figures/localscore/",filename,'_ScoreLocalAll_xi',xi,'.pdf'))
#   par(mfrow=c(4,1), mar=c(5,2,1,1))
#   for (g in chrInfo$chr){
#     plot(fet_df[chr==g,pos],fet_df[chr==g,lindley],xlab=paste("position (chr ",g,")",sep=''),ylab="Score Local", type='l', ylim=c(0,max(fet_df[chr==g,max(lindley)],chrInfo[chr==g, th])))
#     abline(h=chrInfo[chr==g, th], col='grey')
#     abline(v=sigZones[chr==g, beg], col='red', lty=3)
#     abline(v=sigZones[chr==g, end], col='red', lty=3)
#   }
#   dev.off()
# }
# 
# ### In this example there are too many significant regions because the distribution
# #of the p-values is not uniform due to experimental designs


### FST ######

chromfile <- read_tsv("data/R/chrom_meta.txt", show_col_types = F)

alpha = "01"

ls_sig_df <- read_tsv(paste0('./results/localscore/parameterization/pink-chum_minInd0.3_MM3_xi',xi,'_',alpha,'sig.txt')) %>%
  rename(chrName = chr) %>%
  left_join(chromfile,by='chrName') %>%
  select(-chrName)

pink_list <- as.list(Sys.glob("./results/fst/*pink-chum*early-late_minInd0.3*fst.txt"))

# read in all data files that match wildcard
pink_df <- pink_list %>%
  set_names(nm = pink_list) %>%   
  map_dfr(
    ~ read.delim2(.x, header = T, sep = "\t", row.names = NULL, 
                  col.names = c("region", "chrName", "midPos", "Nsites", "Fst"))
  ) %>% 
  select(-c(region, Nsites)) %>%
  mutate(midPos = as.numeric(midPos/10^6),
         Fst = as.numeric(Fst),
         Fst = ifelse(Fst < 0,0,Fst)) %>%
  left_join(chromfile,by='chrName') %>%
  select(-chrName)

fst_localpeaks <- ggplot() +
  geom_point(data = pink_df, aes(x=midPos, y=Fst, color = as_factor(chr))) +
  geom_vline(data = ls_sig_df, aes(xintercept = beg/1e6), color = "firebrick", linetype="dashed",alpha=0.5) +
  geom_vline(data = ls_sig_df, aes(xintercept = end/1e6), color = "blue", linetype="dotted",alpha=0.5) +
  facet_grid(~chr, scales = 'free',space="free") +
  scale_color_manual(values = rep(c("gray50","black"),length(unique(pink_df$chr))))+
  scale_x_continuous(expand = expansion(mult = c(0.003, 0.003)),
                     breaks = seq(20,100,by=20))+
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.001)),
                     breaks = seq(0, 1, by = 0.2),
                     limits = c(0, 1.02)) +
  theme_classic() +
  theme(plot.margin = unit(c(0,0.1,0.1,0.1), "cm"),
        legend.position = "none")

ggsave(paste0("./figures/localscore/pink-chum_fst_localscore_xi",xi,"_sig",alpha,".png"),
       plot = fst_localpeaks, width = 30, height = 8)
