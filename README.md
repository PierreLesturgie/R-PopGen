# Population Genomics with R
----------------------
#### Author: Pierre Lesturgie
----------------------

## Analyse and filter genomics data using R functions
 
This repository provides python scripts to analyse and filter genomic data.
 
Installing the following packages will ensure all functions to work:
- parallel
- vcfR
- ggplot2
- ggthemes
- scales
- ade4
- rlist

## (1) Filtering and extracting info from VCFs
### Missing data
Returns a list with missing data per sample or SNPs or both
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)

    missing.data(GT=NULL,vcf=NULL, SAMPLE=T, SNPs=T)

### Depth
Returns a list with depth of coverage data per sample or SNPs or both
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)

    depth(DP=NULL,vcf=NULL, SAMPLE=T, SNPs=T)
    
### gc content
Returns GC data
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)

    gc.content(GT=NULL, vcf=NULL)

### Heterozygosity
Returns a list with heterozygosity data per sample or SNPs or both
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)

    heterozygosity(GT=NULL,vcf=NULL, SAMPLE=F, SNPs=T)

### Frequency of fixed alleles
Returns a list with frequency of fixed alleles data per sample or SNPs or both
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)

    fixed.alleles <- function(GT=NULL,vcf=NULL, SAMPLE=F, SNPs=T)

### Minor allele frequency
Returns MAF data
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)

    minor.allele.frequency(GT=NULL,vcf=NULL)

### Pick random SNP per locus
Returns a filtered VCF. This is usually meaningful for RAD-like data (i.e., where there are many loci)

    random.snp<-function(vcf)

## (2) Site frequency spectrum (SFS)

### Folded SFS calculation
Returns folded SFS (accross all SNPs or in sliding windows)
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function)
If needing monomorphinc sites (i.e., with a minor allele frequency of 0) use <mono=T>
Possibility of parallelization 

    sfs.folded = function(vcf=NULL,GT=NULL,slide=100000,jump=25000,sliding_window=FALSE,paral=F,mono=F)

### Unfolded SFS calculation
Returns folded SFS (and possibly a derived vcf)
Need to input either a vcf or at GT matrix (obtained from vcfR extract.gt() function) and the sample name of the outgroup

    sfs.unfolded(vcf=NULL, GT=NULL, outgroup, return_derived_vcf=F)

### SFS folding 
Returns folded SFS from unfolded

    fold.SFS(sfs)

### Calculates a correct SFS with missing data
Returns a folded SFS from a vcf and a maximum missing data rate

    corr.sfs.missing.data<-function(vcf,md)

### SFS normalization as in Lapierre et al. (2021)
Returns normalized SFS from sfs (Without monomorphic sites)
Possibility to return the expected norm SFS and expected SFS under and panmictic and constant population scenario

    normalized.expected.SFS(sfs, return.expected.norm=TRUE, return.norm=TRUE, return.expected.SFS=TRUE)

### SFS normalization as in Lapierre et al. (2021) in sliding windows
Returns normalized SFS from sfs (Without monomorphic sites) per windows

    normalized.expected.SFS.sliding.windows(sfs)

### Euclidian distance in sliding windows
Returns Euclidian distance between X and Y per windows (typically two SFS)

    euclidian.distance.sliding.window(X,Y=avgSFS)
    
### Distance SFS in windows to Average SFS 
Input is a data frame with the first three columns being: <CONTIG>, <low_bound_window>, <high_bound_window>
followed by the SFS for each window

    dist.normSFS(data)

### Calculates Watterson's (1979) estimate of genetic diversity
Input is simply a folded SFS (only polymorphic sites)
    watterson<-function(sfs)


### Psi calculation as in Peter & Slatkin (2013, 2015)
Input is a dataframe with number of derived allees per SNP and individual. 
The second script (psi2dist) returns a distance matrix with PSI valyes and significancy

    psi(der_all, bootstrap=1000)
    psi2dist=function(psi,signif=T)


### Genetic distance between two individuals
There are two functions: 

#### (A) : D = 1 - (shared_alleles / 2)
Input is simply a genotype table (obtained using vcfR extract.gt function)

    genetic.distance=function(gt)

#### (B) : Bray-Curtis distance
Input is simply a vcf

    BC.distance.individual(vcf,bootstrap=NULL)

### Plot Stairway Plot with confidence intervals
Input is a **named** list of multiple output dataframe from stairwayplot (*.summary)
Also works with a single dataframe but still needs to be stored in a list. 

    plot.stairway.IC<-function(data_list,var="Ne",output='stairway.pdf',alpha=0.1,cols=NULL,CI=T,xlim=NULL, legend=T,
                           leg_pos=c(0.78, 0.92),ylim=NULL,x.breaks=4,y.breaks=5,ncol.leg=1,by.gen=NULL,
                           vline=NULL,col_vline='grey',ltyp=NULL,ltyp_vline='dashed')

### Neighbouring SFS in ABC sumstat
Computes the N closest SFS to the observed one in one or multiple simulated models 

    neighbours.SFS<-function(sumstat,target,n.closest=100,plot.output=NULL)
    multi.neigh.SFS<-function(sumstat_list,target,n.closest=100,plot.output=NULL)

### Fluctuations in stairway plot output
computes the sum of slopes between discretized intervals from the stairway plot output
    fluctuations.stairway(stairway,N_DISCRETIZE)


# ========================================== #
# >>>> III. Hudson's FST <<<< #
# ========================================== #

gt2alleles = function(gt){
  alleles = unlist(strsplit(gt,"[/|]"))
  alleles[which(alleles==".")] = NA
  while(length(alleles) != length(gt)*2){
    alleles <- c(alleles,NA)
  }
  return(as.numeric(alleles))
}

mpd = function(an,ac){
  n_pairs = an * (an - 1) / 2
  n_same = rowSums(ac * (ac - 1) / 2)
  
  n_diff = n_pairs - n_same
  
  return(n_diff / n_pairs)
}

mpd.between = function(an1,an2,ac1,ac2){
  n_pairs = an1 * an2 
  n_same = rowSums(ac1 * ac2)
  
  n_diff = n_pairs - n_same
  
  return(n_diff / n_pairs)
}

p.FST = function(an1,an2,ac1,ac2){
  if(is.null(dim(ac1))){
    ac1 = t(data.frame(ac1))
    ac2 = t(data.frame(ac2))
  }
  within = (mpd(an1,ac1) + mpd(an2,ac2))/2
  between = mpd.between(an1,an2,ac1,ac2)
  fst_nuc = (between - within) / between
  fst = (sum(between,na.rm = T) - sum(within,na.rm = T))/sum(between,na.rm = T)
  return(list(fst,fst_nuc))
}

sliding.window.fst = function(windows,positions,an1,an2,ac1,ac2){
  temp<-which(positions >= as.numeric(windows[1]) & positions < as.numeric(windows[2]))
  if(length(temp) > 0){
    pairwsie_fst = p.FST(an1 = an1[temp], an2 = an2[temp], 
                         ac1 = ac1[temp,], ac2 = ac2[temp,])
    return(pairwsie_fst[[1]])
  }else{return(NA)}
}

a.counts = function(site,mx){
  c <- c()
  for(i in 1:(mx+1)){
    c <- c(c,length(which(site == (i-1))))
  }
  return(c)
}

### Wrapping function 

fst.hudson <- function(vcf, pop_list,resampling=100,sliding_window=FALSE,slide=NULL,jump=NULL,write=FALSE){
  require(vcfR)
  require(rlist)
  
  GT = extract.gt(vcf,as.numeric = F,convertNA = F)
  
  cat("\n","<< Counting Alleles Per Individual >>", "\n")
  alleles_indv = list()
  pb <- txtProgressBar(min = 0, max = ncol(GT), style = 3)
  for(i in 1:ncol(GT)){
    alleles_indv[[i]] = t(apply(data.frame(GT[,i]),1,gt2alleles))
    setTxtProgressBar(pb, i)
  }
  
  names(alleles_indv) = colnames(GT)
  
  
  alleles_counts = allele_num = list()
  for(i in 1:length(pop_list)){
    
    alleles<-c()
    for(p in 1:length(pop_list[[i]])){
      alleles <- cbind(alleles,alleles_indv[[which(names(alleles_indv)==pop_list[[i]][p])]])
    }
    
    #alleles = t(apply(GT[,pop_list[[i]]],1,gt2alleles))
    alleles <- split(alleles, seq(nrow(alleles)))
    alleles_counts[[i]] = t(mapply(a.counts,alleles,MoreArgs = list(mx = max(as.numeric(unlist(strsplit(max(GT),"[/|]")))))))
    allele_num[[i]] = rowSums(alleles_counts[[i]])
    
  }
  
  if(sliding_window == TRUE){
    positions = getPOS(vcf)
    cat("\n","<< Computing Pairwise FST with Sliding Widows >>", "\n")
  }else{cat("\n","<< Computing Pairwise FST >>", "\n")}
  
  pb <- txtProgressBar(min = 0, max = length(pop_list)*(length(pop_list)-1)/2, style = 3)
  fst_nuc = sliding_win_result = list()
  fst = nms = c()
  r=0
  for(i in 1:(length(pop_list)-1)){
    
    for(j in (i+1):length(pop_list)){
      
      if(sliding_window==TRUE){
        windows<-data.frame(low_bound=seq(1, max(positions) + jump, jump),
                            upper_bound=seq(slide, max(positions) + jump + slide, jump))
        nw.list <- split(windows, seq(nrow(windows)))
        sw = mapply(sliding.window.fst,
                    nw.list,
                    MoreArgs = list(positions = positions,
                                    an1 = allele_num[[i]], 
                                    an2 = allele_num[[j]],
                                    ac1 = alleles_counts[[i]], 
                                    ac2 = alleles_counts[[j]]))
        
        sliding_win_result = list.append(sliding_win_result,
                                         data.frame(windows,
                                                    mid_wind=((as.numeric(windows[,1]) + 
                                                                 as.numeric(windows[,2])) / 2),
                                                    FST=sw))
      }
      
      pairwsie_fst = p.FST(an1 = allele_num[[i]], an2 = allele_num[[j]], 
                           ac1 = alleles_counts[[i]], ac2 = alleles_counts[[j]])
      
      fst = c(fst,pairwsie_fst[[1]])
      fst_nuc = list.append(fst_nuc,pairwsie_fst[[2]])
      nms = c(nms,paste0(names(pop_list)[i],"/",names(pop_list)[j]))
      
      r = r + 1
      setTxtProgressBar(pb, r)
    }
    
  }
  
  names(fst_nuc) = nms
  if(sliding_window == TRUE){names(sliding_win_result) = nms}
  
  fst_obs = data.frame(FST=fst)
  rownames(fst_obs) = nms
  
  if(is.null(resampling)==FALSE){
    cat("\n","<< Computing",resampling,"resampling runs on Pairwise FST >>", "\n")
    
    pb <- txtProgressBar(min = 0, max = resampling, style = 3)
    fst_all<-c()
    for(r in 1:resampling){
      
      fst<-c()
      for(i in 1:(length(pop_list)-1)){
        
        for(j in (i+1):length(pop_list)){
          
          pop_temp <- c(pop_list[[i]],pop_list[[j]])
          s = sample(length(pop_temp))
          
          pop1 = pop_temp[s[1:length(pop_list[[i]])]]
          pop2 = pop_temp[-s[1:length(pop_list[[i]])]]
          
          pop_list_res = list(pop1,pop2)
          
          alleles_counts = allele_num =  list()
          for(f in 1:2){
            
            alleles<-c()
            for(p in 1:length(pop_list_res[[f]])){
              alleles <- cbind(alleles,alleles_indv[[which(names(alleles_indv)==pop_list_res[[f]][p])]])
            }
            
            alleles <- split(alleles, seq(nrow(alleles)))
            alleles_counts[[f]] = t(mapply(a.counts,alleles,MoreArgs = list(mx = max(as.numeric(unlist(strsplit(max(GT),"[/|]")))))))
            #alleles_counts[[f]] = t(apply(alleles,1,a.counts))
            allele_num[[f]] = rowSums(alleles_counts[[f]])
          }
          
          pairwsie_fst = p.FST(an1 = allele_num[[1]], an2 = allele_num[[2]], 
                               ac1 = alleles_counts[[1]], ac2 = alleles_counts[[2]])
          fst = c(fst,pairwsie_fst[[1]])
        }
        
      }
      names(fst) = rownames(fst_obs)
      if(write==T){
        FST=data.frame(fst)
        write.table(t(FST),paste0("run_",r,".fst"),row.names = F)
      }
      
      
      fst_all <- rbind(fst_all,fst)
      
      setTxtProgressBar(pb, r)
    }
    
  }else{fst_all = NULL}
  
  
  result = list(fst_obs,fst_nuc,fst_all,sliding_win_result)
  names(result) = c("pairwise_fst","pairwise_fst_nuc","pairwise_fst_resampling","pairwise_fst_sliding_window")
  
  return(result)
}

#create a distance matrix from geo coordinates (lon, lat)
distance.geo<-function(map, unit="km"){
  require(geosphere)
  
  mat<-matrix(NA, nrow = length((map$x)), ncol=length((map$x)))
  rownames(mat)<-rownames(map)
  colnames(mat)<-rownames(map)
  for (i in 1:length((map$x)-1)){
    for (j in i:length((map$x))){
      if (j==i){
        mat[i,j]=0
      }
      else
      {
        dgeo<-distGeo(map[c(i,j),])
        mat[i,j]<-dgeo[1]
        mat[j,i]<-dgeo[1]
      }
      
    }
    
  }
  
  if (unit=="m"){
    mat_m<-as.dist(mat)
    print("distances in meters")
    return(mat_m)
  }
  if (unit=="km"){
    mat_km<-mat/1000
    mat_km<-as.dist(mat_km)
    
    print("distances in kilometers")
    return(mat_km)
  }
  
}

dist.fst<-function(df_fst, popmap=NULL, signif=F, pv=0.01, sep_pop="/"){
  if (is.null(popmap)){
    popnames<-rownames(df_fst)
    popnames<-as.character(popnames)
    poplist<-unlist(strsplit(popnames, split = sep_pop))
    dbl<-which(duplicated(poplist))
    pop<-poplist[-dbl]
  } else {
    if (is.null(dim(popmap))==F){
      pop<-levels(popmap[,2])
    }else {pop<-popmap}
    
  }
  
  npop<-length(pop)
  
  matrix_fst<-matrix(NA, ncol=npop,nrow=npop)
  rownames(matrix_fst)<-pop
  colnames(matrix_fst)<-pop
  
  #Inset du df initial
  g<-c()
  g[1]<-1
  for (p in 2:(npop-1)){
    g[p]<-g[p-1]+npop-(p-1)
  }
  
  #Insert a partir de ligne 2
  npop
  d<-c(2:npop)
  
  for (i in 1:(npop-1)){
    h<-g[i]
    for (j in d[i]:npop){
      q<-j-d[i]
      matrix_fst[j,i]<-df_fst[h+q,1] ######ATTENTION FIRST COL : FST; 2ND : 2.5%
      if(signif==T){
        if (df_fst[h+q,1]>df_fst[h+q,4]){matrix_fst[i,j]<-paste0("p<",as.character(pv))}else {matrix_fst[i,j]<-"NS"}
      }
    }
  }
  if (signif==T){
    for (i in 1:length(matrix_fst[,1])){
      matrix_fst[i,which(is.na(matrix_fst[i,]))]<-0
    }
    matrix_fst<-as.data.frame(matrix_fst)
  }else {matrix_fst<-as.dist(matrix_fst)}
  return(matrix_fst)
}

ibd.fst<-function(dgeo,dfst){
  require(ade4)
  dfst<-dfst/(1-dfst)
  ibd <- mantel.randtest(dfst,dgeo, nrepet = 999)
  #print(ibd)
  
  plot(dgeo, dfst, pch=20,cex=.5, ylab="FST/(1-FST)", xlab="Geographic distance (km)")
  d<-lm(dfst~dgeo)
  abline(d$coefficients)
  title("Isolation by distance")
  
  fin<-list(ibd, dgeo, dfst)
  names(fin)<-c("IDB","geo_dist","FST_dist_inv")
  return(fin)
}

plot.fst<-function(FST,output="",gradmid=F){
  require(qgraph)
  require(reshape2)
  require(ggplot2)
  require(ggthemes)
  
  if(class(FST)=="dist"){matrix_fst<-as.matrix(FST)}else if (class(FST)=="matrix"){matrix_fst<-FST}
  
  
  #### THE FOLLOWING LOOP CHANGES NEGATIVE VALUES IN POSITIVE ONES (ARTEFACT OF METHOD USED TO COMPUTE FST)
  for (i in 1:length(matrix_fst[,1])){
    n<-which(matrix_fst[i,]<0)
    matrix_fst[i,n]<-matrix_fst[i,n]*(-1)
  }
  
  similarity_matrix <- 1/matrix_fst # one over, as qgraph takes similarity matrices as input
  
  #this is a network of similarity
  pdf(paste0(output,'similarity_network.pdf'))
  qgraph(similarity_matrix, layout='groups')
  dev.off()
  
  
  #This is now to display a graph with pairwise-fst values
  melted_dfst <- melt(matrix_fst)
  #head(melted_dfst)
  
  for(i in 1:nrow(melted_dfst)){
    if(melted_dfst$value[i]==0){melted_dfst$value[i]<-NA}
  }
  
  # col "royalblue4" and "lightcyan" are good
  g1<-ggplot(data = melted_dfst, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() + labs(x="", y="", fill=expression(F[ST])) +
    scale_fill_gradient(high="#D53E4F", low="#FFFFBF",na.value = "white") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.text.y = element_text(),
          text = element_text(size = 40)) + theme_few()
  if(gradmid==T){
    g1<-ggplot(data = melted_dfst, aes(x=Var1, y=Var2, fill=value)) +
      geom_tile() + labs(x="", y="", fill=expression(F[ST])) +
      scale_fill_gradient2(high="black",mid="#D53E4F", low="#FFFFBF",na.value = "white") +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.text.y = element_text(),
            text = element_text(size = 40)) + theme_few()
  }
  
  
  
  
  print(g1)
  ggsave(filename = paste0(output,"pairwise_FST.pdf"),g1,height = 5,width = 6)
  return(g1)
}


# ========================================== #
# >>>> IV. OTHER USEFUL FUNCTIONS <<<< #
# ========================================== #

euclidian.distance = function(X,Y){
  if(length(X)!=length(Y)){
    stop("X and Y don't have the same length")
  }
  Z=c()
  for(i in 1:length(X)){
    Z = c(Z,(X[i]-Y[i])^2)
  }
  Z = sqrt(sum(Z,na.rm = T))
  return(Z)
}

index.n.inter <- function(fst,phylogenetic_scale="species",threshold=0){
  
  fst[,paste0(phylogenetic_scale)] <- as.character(fst[,paste0(phylogenetic_scale)])
  
  nais = which(is.na(fst[,paste0(phylogenetic_scale)]))
  if(length(nais > 0)){
    fst = fst[-nais,]
  }
  
  if(length(levels(as.factor(fst[,paste0(phylogenetic_scale)])))!=nrow(fst)){
    cat("","Summing by prey",phylogenetic_scale, "\n")
    final <-c()
    pb <- txtProgressBar(min = 0, max = length(levels(as.factor(fst[,paste0(phylogenetic_scale)]))), style = 3)
    for(i in 1:length(levels(as.factor(fst[,paste0(phylogenetic_scale)])))){
      temp = fst[which(fst[,paste0(phylogenetic_scale)]==levels(as.factor(fst[,paste0(phylogenetic_scale)]))[i]),]
      temp = temp[,-c(1:10)]
      temp = colSums(temp)
      final = rbind(final,temp)
      setTxtProgressBar(pb, i)
    }
    
    interactions = t(final)
  }else{
    interactions = t(fst[,-c(1:10)])
  }
  rnmes = rownames(interactions)[-nrow(interactions)]
  
  d = matrix(unlist(strsplit(rnmes,split = "_")),ncol = 5,byrow = T)[,c(2,3)]
  d = paste0(d[,1],"_",d[,2])
  
  rownames(interactions) = c(d,"SUM")
  
  cat("\n","Summing by consumer species", "\n")
  interactions_2<-c()
  pb <- txtProgressBar(min = 0, max = length(levels(as.factor(d))), style = 3)
  for(i in 1:length(levels(as.factor(d)))){
    temp = levels(as.factor(d))[i]
    #print(temp)
    if(length(which(rownames(interactions)==temp))>1){
      temp2 = colSums(interactions[which(rownames(interactions)==temp),])
      interactions_2 = rbind(interactions_2,c(temp,temp2))
    }else{
      interactions_2 = rbind(interactions_2,c(temp,interactions[which(rownames(interactions)==temp),]))
    }
    setTxtProgressBar(pb, i)
  }
  
  interactions_2 <- rbind(interactions_2, c("SUM",interactions[nrow(interactions),]))
  
  cat("\n","Counting interactions", paste0("(threshold = ",threshold,")"),"\n")
  interactions_clean<-c()
  pb <- txtProgressBar(min = 2, max = ncol(interactions_2), style = 3)
  for(i in 2:ncol(interactions_2)){
    temp = interactions_2[,i]
    summ = temp[nrow(interactions_2)]
    if(as.numeric(summ)>threshold){
      temp2 = as.numeric(temp[-(nrow(interactions_2))]) / as.numeric(summ)
    }else{
      temp2 = as.numeric(temp[-(nrow(interactions_2))])
    }
    interactions_clean <- cbind(interactions_clean,c(temp2))
    setTxtProgressBar(pb, i)
  }
  
  number_of_interactions=c()
  for(i in 1:nrow(interactions_clean)){
    temp = interactions_clean[i,]
    number_of_interactions = c(number_of_interactions,length(which(temp>0)))
  }
  
  interactions_clean <-data.frame(species=interactions_2[-nrow(interactions_2),1],N_inter=number_of_interactions)
  return(interactions_clean)
}
