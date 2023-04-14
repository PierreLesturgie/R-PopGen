#################################################################
#---------------------------------------------------------------#
#----------- FUNCTIONS RELATED TO GENETIC DIVERSITY ------------#
#---------------------------------------------------------------#
#################################################################


# 1. Computes Folded Site Frequency Spectrum from VCF
# >>> There is a parallel option that might not work on windows machine
# >>> Sliding window option requires having data from a single chromosome so far
# >>> If mono=T, first class is number of monomorphic sites (useful for "complete" VCFs)
sfs.folded = function(vcf=NULL,GT=NULL,slide=100000,jump=25000,sliding_window=FALSE,paral=F,mono=F){
  
  if(is.null(GT)){
    require(vcfR)
    
    GT<-extract.gt(vcf, element = "GT", mask = FALSE,
                   as.numeric=F,return.alleles = FALSE, 
                   IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  }
  
  nind<-ncol(GT)
  
  minor.allele.count<-function(gt){
    alleles <- unlist(strsplit(gt,split = "[/|]"))
    A <- length(which(alleles=="0"))
    B <- length(which(alleles=="1"))
    return(min(A,B))
  }
  
  GT.list<-split(GT, seq(nrow(GT)))
  
  if(paral==T){
    require(parallel)
    cores<-(detectCores()-1)
    cl <- makeCluster(mc <- getOption("mc.cores", cores))
    q<-mcmapply(minor.allele.count,GT.list,mc.cores = c(detectCores()-1))
    stopCluster(cl)
  }else{q<-mapply(minor.allele.count,GT.list)}
  
  
  
  sfs<-matrix(NA,nrow = 1, ncol=(nind+1))
  for (j in 0:nind){
    sfs[1,j+1]<-length(which(q==j))
  }
  
  if(mono==F){
    sfs = sfs[-1]
  }
  
  if(sliding_window==T){
    
    if(is.null(vcf)==F){
      positions<-getPOS(vcf)
    }else{
      positions<-as.numeric(matrix(unlist(strsplit(rownames(GT),"_")),byrow = T,ncol=3)[,3])
    }
    
    sliding.sfs<-function(wind,nind,positions){
      sfs_temp<-matrix(NA,nrow = 1, ncol=nind)
      temp<-which(positions >= as.numeric(wind[1]) & positions < as.numeric(wind[2]))
      if(length(temp)>0){
        q_temp<-q[temp]
        for (j in 1:nind){
          sfs_temp[1,j]<-length(which(q_temp==j))
        }
      }
      return(sfs_temp)
    }
    
    nw<-data.frame(low_bound=seq(1, max(positions) + jump, jump),
                   upper_bound=seq(slide, max(positions)+jump+slide, jump))
    
    nw.list <- split(nw, seq(nrow(nw)))
    if(paral==T){
      cl <- makeCluster(mc <- getOption("mc.cores", cores))
      sfs_win<-t(mcmapply(sliding.sfs,nw.list,MoreArgs = list(nind,positions),mc.cores = c(detectCores()-1)))
      stopCluster(cl)
    }else{sfs_win<-t(mapply(sliding.sfs,nw.list,MoreArgs = list(nind,positions)))}
    
    rm(nw.list)
    rownames(sfs_win)<-(nw[,2]+nw[,1])/2
  }
  
  if(sliding_window==T){
    final<-list(sfs,sfs_win) ; names(final) <- c("sfs_global","sliding_sfs")
    return(final)
  }else{return(sfs)}
  
}

# 2. Function to compute an unbiased SFS in case of missing data 
# >>> md is the missing data rate allowed
sfs.folded.missing.data<-function(vcf,md){
  
  require(vcfR)
  
  haploid<-ncol(vcf@gt)-1
  diploid<-haploid*2
  
  haps<-t(extract.haps(vcf,unphased_as_NA = F))
  
  sequences<-matrix(1,nrow = dim(haps)[1],ncol = dim(haps)[2])
  for(i in 1:ncol(haps)){
    na<-which(is.na(haps[,i]))
    t<-table(haps[,i])/sum(table(haps[,i]))
    if(as.numeric(t[1])==0.5){
      maj<-names(t)[1]
    }else{maj<-names(t)[(which(t==max(t)))]}
    sequences[which(haps[,i]==maj),i]<-0
    sequences[which(is.na(haps[,i])),i]<-NA
  }
  
  
  
  sequences_new<-matrix(ncol = ncol(sequences),nrow = c(diploid-md*diploid))
  for(i in 1:ncol(sequences)){
    n<-which(is.na(sequences[,i]))
    p<-(length(n)/nrow(sequences))
    if(p!=md){
      if(p==0){
        d<-sample((diploid-length(n)),(md)*diploid,replace = F)
        d<-d[order(d)]
        sequences_new[,i]<-sequences[-d,i]
      }else{
        sequence_temp<-sequences[-n,i]
        d<-sample((diploid-length(n)),((md)*diploid-p*diploid),replace = F)
        d<-d[order(d)]
        sequences_new[,i]<-sequence_temp[-d]
      }
    }else{
      sequences_new[,i]<-sequences[-n,i]
    }
  }
  
  
  SNPS<-c()
  for(i in 1:ncol(sequences_new)){
    s<-data.frame(table(sequences_new[,i]))
    k1<-s[which(s$Var1==1),2];if(length(k1)==0){k1<-0}
    k0<-s[which(s$Var1==0),2];if(length(k0)==0){k0<-0}
    if(k0 < k1){
      k<-s[which(s$Var1==0),2]
    }else{
      k<-s[which(s$Var1==1),2]
    }
    if(length(k)>0){
      SNPS<-c(SNPS,k)
    }else{SNPS<-c(SNPS,0)}
  }
  new_sfs<-as.numeric(table(SNPS))[-1]
  
  return(new_sfs)
}

# 3. Computes Unfolded Site Frequency Spectrum from VCF with an outgroup
sfs.unfolded<-function(vcf=NULL, GT=NULL, outgroup, return_derived_vcf=F){
  
  if(is.null(GT)){
    require(vcfR)
    GT<-extract.gt(vcf, element = "GT", mask = FALSE,
                   as.numeric=F,return.alleles = FALSE, 
                   IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  }
  
  outgroup_idx<-which(colnames(GT)==as.character(outgroup))
  
  # Remove heterozygous sites
  idx <- as.numeric(which(GT[,outgroup_idx]== "0/1" | GT[,outgroup_idx]== "0|1"))
  
  if (length(idx)>0) {
    GT<-GT[-idx,]
  }
  
  # function for computing derived frequences related to the outgroup
  derived.allele.frequency<-function(gt){
    temp<-gt[-outgroup_idx]
    N<-length(temp)*2
    freq <- N - ((length(which(temp==gt[outgroup_idx])))*2 + (length(which(temp=="0/1" | temp=="0|1"))))
    return(freq)
  }
  
  SNPs_freq<-apply(GT,MARGIN = 1,derived.allele.frequency)
  
  sfs=c()
  for (i in 0:((ncol(GT)-1)*2)){
    n<-length(which(SNPs_freq==i))
    sfs<-c(sfs,n)
  }
  
  sfs<-matrix(sfs,nrow = 1) ; colnames(sfs) <- c(0:((ncol(GT)-1)*2))
  
  if(return_derived_vcf==T){
    return(list(sfs,vcf[-idx,-c(outgroup_idx+1)]))
  }else{return(sfs)}
}

# 4. Folds an unfolded SFS
fold.SFS<-function(sfs){
  d<-c()
  for(i in 1:(((length(sfs)+1)/2))-1){
    d = c(d, sfs[i] + sfs[length(sfs) + 1 - i])
  }
  d<-c(d, sfs[((length(sfs) + 1) / 2)])
  return(d)
}

# 5. Computes the normalized SFS as in Lapierre et. al 2017 from SFS
# >>> can also return the expected nSFS and SFS under a panmictic constant demographic model. 
normalized.expected.SFS <- function(sfs, return.expected.norm=FALSE, return.norm=TRUE, return.expected.SFS=FALSE){
  require(rlist)
  
  RES = list()
  norma_sfs = sfs / sum(sfs)
  ind = length(norma_sfs)
  
  n_sfs <- c()
  for(i in 1:(ind-1)){
    n_sfs<-c(n_sfs,as.numeric(norma_sfs[i])* i * ((2*ind-i))/(2*ind))
  }
  n_sfs <- c(n_sfs,as.numeric(norma_sfs[ind]*ind))
  
  X = seq(1,ind,1)/(2*ind)
  
  
  if(return.norm==TRUE){
    temp = data.frame(X=X,n_sfs=n_sfs)
    RES = list.append(RES,temp)
    names(RES)[length(RES)] = "Norm_SFS"
  }
  
  if(return.expected.norm==TRUE | return.expected.SFS==TRUE){
    expected_norm = rep(sum(n_sfs) / ind, ind)
    
    if(return.expected.norm==TRUE){
      temp = data.frame(X=X,n_sfs=expected_norm)
      RES = list.append(RES,temp)
      names(RES)[length(RES)] = "Expected_Norm_SFS"
    }
    
    if(return.expected.SFS==TRUE){
      expected_norma <- c()
      for(i in 1:(ind-1)){
        expected_norma = c(expected_norma,expected_norm [i] / (i * ((2*ind-i))/(2*ind)))
      }
      expected_norma <- c(expected_norma,as.numeric(expected_norm [ind]/ind))
      expected = sum(sfs) * expected_norma
      RES = list.append(RES,expected)
      names(RES)[length(RES)] = "Expected_SFS"
    }
  }
  return(RES)
}

# 6. Computes Psi as in Peter & Slatkin (), the directionality index from a 3D SFS with three individuals
# >>> One individual is the outgroup, the other two are from two populations
# >>> Computes a resampling of alleles to assess significancy of Psi
psi.from.3D.SFS<-function(SFS=NULL,pop,directory,bootstrap=999){
  
  npop<-length(pop)
  
  pair_populations<-c()
  for(i in 1:(npop-1)){
    pair_populations<-c(pair_populations,paste0(rep(pop[i],npop-i),"_",pop[(i+1):npop]))
  }
  
  
  matrice_psi<-data.frame(PSI=matrix(NA,nrow = npop,ncol=npop),row.names = pop);colnames(matrice_psi)<-pop
  matrice_signif<-data.frame(PSI=matrix(NA,nrow = npop,ncol=npop),row.names = pop);colnames(matrice_signif)<-pop
  
  nsites_der<-list()
  for(i in 1:length(pair_populations)){
    if(is.null(SFS)==T){
      temp<-read.table(paste0(directory,pair_populations[i],"_3dsfs.sfs"))
    }else{temp<-SFS[[i]]}
    
    # temp[6] is pop a hetero for fixed (homo) pop b >> the derived alleles for pop a
    # exactly it is the entry (0,1,2) of the 3d sfs : derived SNP fixed in pop3 + hetero in pop2 (pop1 is outgroup)
    # temp[8] is pop b hetero for fixed (homo) pop a >> the derived alleles for pop b
    # exactly it is the entry (0,2,1) of the 3d sfs : derived SNP fixed in pop2 + hetero in pop3 (pop1 is outgroup)
    # temp[1] is the total number of sites
    psi_obs<-(round(temp[8]) - round(temp[6]))/temp[1]
    
    nsites<-round(temp[8]) + round(temp[6])
    nsites_der[[i]]<-cbind(temp[8],temp[6],temp[1],psi_obs) 
    colnames(nsites_der[[i]])<-c(unlist(strsplit(pair_populations[i],"_")),"n_sites","psi")
    names(nsites_der)[i]<-pair_populations[i]
    
    psi_distrib<-c()
    for(j in 1:bootstrap){
      test_temp1<-rbinom(n = 1,size = as.numeric(nsites),prob = 0.5)
      test_temp2<-as.numeric(nsites) - test_temp1
      psi_temp<-(test_temp1 - test_temp2)/temp[1]
      psi_distrib<-c(psi_distrib,as.numeric(psi_temp))
    }
    
    if(psi_obs<0){
      h<-which(psi_distrib<=as.numeric(psi_obs))
    }else{
      h<-which(psi_distrib>=as.numeric(psi_obs))
    }
    
    
    p=length(h)/(bootstrap+1)
    
    # so if psi < 0 : a has less derived alleles than b. 
    # so pop b 
    ppp<-c(unlist(strsplit(pair_populations[i],split = "_")))
    a<-which(rownames(matrice_psi)==ppp[1])
    b<-which(rownames(matrice_psi)==ppp[2])
    matrice_psi[a,b]<-psi_obs;matrice_psi[b,a]<-psi_obs*(-1)
    matrice_signif[a,b]<-p;matrice_signif[b,a]<-p
  }
  
  for (i in 1:length(matrice_psi[,1])){
    matrice_psi[i,which(is.na(matrice_psi[i,]))]<-0
  }
  list_fin<-list(matrice_psi,matrice_signif,nsites_der);names(list_fin)<-c("Psi","signif","nsites_der")
  write.table(matrice_psi,paste0(directory,"/psi.txt"))
  write.table(matrice_signif,paste0(directory,"/psi_signif.txt"))
  return(list_fin)
}

# 7. Function to plot the IICR through time estimated by the stairwayplot
# >>> in data_list one can put many different df of stairwayplots to plot
# >>>  the list must always be named
plot.stairway.IC<-function(data_list,var="Ne",output='stairway.pdf',alpha=0.1,cols=NULL,CI=T,xlim=NULL, legend=T,
                           leg_pos=c(0.78, 0.92),ylim=NULL,x.breaks=4,y.breaks=5,ncol.leg=1,by.gen=NULL,
                           vline=NULL,col_vline='grey',ltyp=NULL,ltyp_vline='dashed'){
  
  require(ggplot2)
  require(ggthemes)
  
  if(var=="theta"){
    h<-"theta"
    index<-c(6,3:5)
  }else if(var=="Ne75"){
    index<-c(6,7,10:11)
    h<-"Ne"
  }else{
    index<-c(6:9)
    h<-"Ne"
  }
  
  if(length(data_list)==1){
    all<-data_list[[1]]
    all<-all[,c(index)]
    single_sp=T
    all<-cbind(all,rep(names(data_list),nrow(all)))
  }else{
    single_sp=F
    all<-c()
    nn<-c()
    for(i in 1:length(data_list)){
      all<-rbind(all,data_list[[i]][,c(index)])
      nn<-c(nn,rep(names(data_list)[i],nrow(data_list[[i]])))
    }
    all<-data.frame(cbind(all,nn))
  }
  
  colnames(all)<-c("Time","temp","lower","upper","Species")
  
  all$Species<-factor(all$Species,levels=names(data_list))
  
  if(is.null(by.gen)==F){
    all[,1]=all[,1]/by.gen
    h2='Generations Before Present'
  }else{h2='Years Before Present'}
  
  if(is.null(ltyp)==F){
    g<-ggplot(data = all, aes(x=Time,y=temp, col=Species)) + geom_line(aes(linetype=Species)) + scale_linetype_manual(values = ltyp)
  }else{
    g<-ggplot(data = all, aes(x=Time,y=temp, col=Species)) + geom_line()
  }
  
  g<-g + ggtitle("") + ylab(paste0(h)) + xlab(paste0(h2)) + theme_few() +
    theme(axis.title.y= element_text(angle = 0,vjust = 0.5)) + 
    scale_x_continuous(n.breaks = x.breaks, labels = comma) + 
    scale_y_continuous(n.breaks = y.breaks,labels = comma)
  if(is.null(cols)==F){
    g<-g + scale_color_manual(values=cols) +
      scale_fill_manual(values=cols)
  }
  
  
  if(CI==T){
    
    g<-g + geom_ribbon(aes(ymin=lower,ymax=upper,
                           col=NULL,fill=Species),alpha=alpha)
  }
  if(legend==FALSE){g <- g + theme(legend.position = "none")}else{
    g <- g + theme(legend.position = leg_pos,
                   legend.title = element_blank()) + guides(col = guide_legend(ncol = ncol.leg))
  }
  if(is.null(vline)==F){
    
    daa<-data.frame(c(vline),c(as.character(vline)))
    colnames(daa)<-c("Time","Grp")
    daa$Grp<-as.factor(daa$Grp)
    g<-g+geom_vline(data = daa,aes(xintercept = Time,colour = Grp),show.legend = F,
                    color=c(col_vline),linetype=c(ltyp_vline))
  }
  
  if(is.null(output)==F){
    ggsave(paste0(output),g)
  }
  return(g)
  
}

# 8. Function to plot the normalized sfs from SFS
# >>> in data_list one can put many different sfs to plot
# >>>  the list must always be named
plot.spettro.multple.mods<-function(sfs_list,leg.pos=c(0.8,0.2),cols=cols,ylim=NULL,ncol.leg=2){
  
  require(ggplot2)
  require(ggthemes)
  
  final<-c()
  for(i in 1:length(sfs_list)){
    SFS_temp<-sfs_list[[i]]
    norm_SFS_temp<-normalized.expected.SFS(SFS_temp)$Norm_SFS
    norm_SFS_temp<-cbind(norm_SFS_temp,rep(names(sfs_list)[i],nrow(norm_SFS_temp)))
    colnames(norm_SFS_temp)<-c("asse_x","eta_2","temp")
    final<-rbind(final,norm_SFS_temp)
  }
  
  
  final<-as.data.frame(final)
  final$asse_x<-as.numeric(as.character(final$asse_x))
  final$eta_2<-as.numeric(as.character(final$eta_2))
  
  gg<-ggplot(data=final,aes(x=asse_x,y=eta_2,col=temp)) + geom_line() +
    ggtitle("") + theme_few() +
    theme(axis.title.y= element_text(angle = 0,vjust = 0.5),
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),legend.title = element_blank()) + 
    scale_discrete_manual(values = c(cols),aesthetics = "col") +
    theme(legend.position = leg.pos,legend.box = "horizontal") + 
    scale_x_continuous(n.breaks = 4, labels = comma)+
    scale_y_continuous(n.breaks = 4, labels = comma) + 
    ylab(expression(SFS[NORM])) + xlab('i/2N') + guides(col = guide_legend(ncol = ncol.leg))
  if(is.null(ylim)==F){
    gg<-gg + ylim(c(ylim))
  }
  return(gg)
}


# 9. Functions to plot the closest simulated SFS to the observed SFS
# >>> implemented also if multiple dataset (multi)
neighbours.SFS<-function(sumstat,target,n.closest=100,plot.output=NULL){
  
  require(ggplot2)
  require(ggthemes)
  
  distance<-c()
  for(i in 1:nrow(sumstat)){
    tmp<-as.data.frame(rbind(target,sumstat[i,]))
    distance<-c(distance,dist(tmp))
  }
  
  neighbours<-c()
  data_temp<-distance
  while(length(neighbours) < n.closest){
    tmp1<-which(data_temp==min(data_temp))
    tmp2<-data_temp[tmp1]
    neighbours<-c(neighbours,tmp2)
    data_temp<-data_temp[-tmp1]
  }
  
  neigh_ss<-c()
  for(i in 1:length(neighbours)){
    tmp<-which(distance==neighbours[i])
    neigh_ss<-rbind(neigh_ss,sumstat[tmp,])
  }
  
  p = priors[as.numeric(rownames(neigh_ss)),]
  
  hist(p[,1])
  
  normsfs<-c()
  for(i in 1:nrow(neigh_ss)){
    normsfs<-cbind(normsfs,calcul_normalized_foldedSFS(neigh_ss[i,])[,2])
  }
  
  low.b<-c()
  high.b<-c()
  for(i in 1:nrow(normsfs)){
    low.b<-c(low.b,min(normsfs[i,]))
    high.b<-c(high.b,max(normsfs[i,]))
  }
  
  final<-as.data.frame(cbind(calcul_normalized_foldedSFS(target),low.b,high.b))
  colnames(final)<-c("eta_2","Norm_SFS","low.b","high.b")
  
  
  g<-ggplot(data = final, aes(x=eta_2,y=Norm_SFS)) + geom_line()+ 
    ggtitle("") + ylab(expression(SFS[NORM])) + xlab('i/2N') + theme_few() +
    scale_discrete_manual(aesthetics = "col",values="black") + 
    theme(axis.title.y= element_text(angle = 0,vjust = 0.5),
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + 
    geom_ribbon(aes(ymin=low.b,ymax=high.b,col=NULL),alpha=0.2)
  print(g)
  
  if(is.null(plot.output)==F){
    ggsave(paste0(plot.output),g)
  }
  return(final)
}

multi.neigh.SFS<-function(sumstat_list,target,n.closest=100,plot.output=NULL){
  
  require(ggplot2)
  require(ggthemes)
  
  neigh<-list()
  for(i in 1:length(sumstat_list)){
    neigh[[i]]<-neighbours.SFS(sumstat = sumstat_list[[i]],target=target,n.closest=n.closest,plot.output=NULL)
  }
  
  all_models<-c()
  for(i in 1:length(neigh)){
    temp<-neigh[[i]]
    temp<-cbind(temp,c(rep(names(sumstat_list)[i],length(temp[,1]))));colnames(temp)[5]<-"Model"
    all_models<-rbind(all_models,temp)
  }
  
  
  
  
  g<-ggplot(data = all_models, aes(x=eta_2,y=Norm_SFS)) + geom_line()+
    ggtitle("") + ylab(expression(SFS[NORM])) + xlab('i/2N') + theme_few()  +
    theme(axis.title.y= element_text(angle = 0,vjust = 0.5),
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), legend.position = c(0.8,0.8), legend.title = element_blank()) +
    geom_ribbon(aes(ymin=low.b,ymax=high.b,col=Model,fill=Model),alpha=0.15,data=all_models)
  print(g)
  
  if(is.null(plot.output)==F){
    ggsave(paste0(plot.output),g)
  }
  return(all_models)
  
}

# 10. Function to compute watterson's () estimator of Theta from folded SFS
# >>> implemented also if multiple dataset (multi)
watterson<-function(sfs){
  sample_size = length(sfs)*2
  h=0
  for (j in 1:(sample_size-1)){
    h=h+1/j
  }
  return(sum(as.numeric(sfs))/h) ####SEE FORMULA 
}