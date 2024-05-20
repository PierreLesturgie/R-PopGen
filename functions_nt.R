#################################################################
#---------------------------------------------------------------#
#------------------ BioInfo & PopGen scripts -------------------#
#---------------------------------------------------------------#
#################################################################

## ****** author  = Pierre Lesturgie ****** ##
## ****** version = 2.0 ****** ##
## ****** email = pierrelesturgie@outlook.fr ****** ##
## ****** date = 13.03.2023 ****** ##

require(parallel)
require(vcfR)
require(ggplot2)
require(ggthemes)
require(scales)
require(ade4)
require(rlist)




# ===================================== #
# >>>> I. Site frequency spectrum <<<<< #
# ===================================== #

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

fold.SFS<-function(sfs){
  d<-c()
  for(i in 1:(((length(sfs)+1)/2))-1){
    d = c(d, sfs[i] + sfs[length(sfs) + 1 - i])
  }
  d<-c(d, sfs[((length(sfs) + 1) / 2)])
  return(d)
}

corr.sfs.missing.data<-function(vcf,md){
  
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

normalized.expected.SFS <- function(sfs, return.expected.norm=TRUE, return.norm=TRUE, return.expected.SFS=TRUE){
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

normalized.expected.SFS.sliding.windows = function(sfs){
  
  norma_sfs = sfs / sum(sfs)
  ind = length(norma_sfs)
  
  n_sfs <- c()
  for(i in 1:(ind-1)){
    n_sfs<-c(n_sfs,as.numeric(norma_sfs[i])* i * ((2*ind-i))/(2*ind))
  }
  n_sfs <- c(n_sfs,as.numeric(norma_sfs[ind]*ind))
  
  X = seq(1,ind,1)/(2*ind)
  
  
  return(n_sfs)
}

euclidian.distance.sliding.window = function(X,Y=avgSFS){
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

dist.normSFS = function(data){
  
  nsfs<-t(apply(data[,-c(1:3)],MARGIN = 1,FUN = normalized.expected.SFS.sliding.windows))
  avgSFS=colMeans(nsfs,na.rm = T)
  
  di=apply(nsfs,MARGIN = 1,FUN = euclidian.distance.sliding.window)
  res = data.frame(CHR=data$V1,POS=data$V2,distSFS=di)
  
  return(res)
  
}

watterson<-function(sfs){
  sample_size = length(sfs)*2
  h=0
  for (j in 1:(sample_size-1)){
    h=h+1/j
  }
  return(sum(as.numeric(sfs))/h) ####SEE FORMULA 
}

fix.ij = function(d){
  df_res = c()
  d = d[-c(1:2)]
  pb <- txtProgressBar(min = 0, max = ncol(d), style = 3)
  for(i in 1:ncol(d)){
    for(j in i:ncol(d)){
      if(i!=j){
        d[,c(i,j)]
        
        f12 = length(which(d[,i]==2 & d[,j]==1))
        f21 = length(which(d[,i]==1 & d[,j]==2))
        f11 = length(which(d[,i]==1 & d[,j]==1))
        f22 = length(which(d[,i]==2 & d[,j]==2))
        SUM=sum(f11,f12,f21,f22)
        df_temp = data.frame(Indi=colnames(d)[i],Indj=colnames(d)[j],f12=f12,f21=f21,SUM,
                             PSI12=((f21-f12)/SUM))
        df_res=rbind(df_res,df_temp)
      }
    }
    setTxtProgressBar(pb, i)
  }
  return(df_res)
}

psi = function(der_all, bootstrap=1000){
  psi_temp = fix.ij(der_all)
  
  if(is.null(bootstrap)==FALSE){
    psi=c()
    pb <- txtProgressBar(min = 1, max = nrow(psi_temp), style = 3)
    for(i in 1:nrow(psi_temp)){
      temp = psi_temp[i,]
      psi_distrib<-c()
      for(j in 1:bootstrap){
        test_temp1<-rbinom(n = 1,size = temp$f12+temp$f21,prob = 0.5)
        test_temp2<-temp$f12 + temp$f21 - test_temp1
        psi_boot<-(test_temp1 - test_temp2)/temp$SUM
        psi_distrib<-c(psi_distrib,as.numeric(psi_boot))
      }
      qt=quantile(psi_distrib, prob=c(0.005, 0.995))
      if(temp$PSI12 < qt[2] & temp$PSI12 > qt[1]){
        p = "NS"
      }else{
        p = "p <= 0.01"
      }
      res = data.frame(temp,pval=p)
      psi=rbind(psi,res)
      setTxtProgressBar(pb, i)
    }
  } else{
    psi = psi_temp
  }
  return(psi)
  
  
  
}

psi2dist=function(psi,signif=T){
  samples=c(psi$Indi,psi$Indj)
  samples=samples[-which(duplicated(samples))]
  matrice=matrix(NA,nrow = length(samples), ncol = length(samples))
  rownames(matrice) = colnames(matrice) = samples
  for(i in 1:length(samples)){
    for(j in i:length(samples)){
      if(i!=j){
        matrice[i,j] = psi$PSI12[which(psi$Indi==samples[i] & psi$Indj==samples[j])]
        if(signif==T){
          matrice[j,i] = psi$pval[which(psi$Indi==samples[i] & psi$Indj==samples[j])]
        }else{
          matrice[j,i] = - psi$PSI12[which(psi$Indi==samples[i] & psi$Indj==samples[j])]
        }
      }
    }
  }
  return(matrice)
}

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

fun_geno_mod<-function(data){
  for (i in 1:length(data)) {
    if (data[i]=="./." | data[i]==".|." | is.na(data[i])==T) {data[i]<-NaN}
    if (data[i]=="0/0" | data[i]=="0|0") {data[i]<-0}
    if (data[i]=="0/1" | data[i]=="0|1") {data[i]<-1}
    if (data[i]=="1/1" | data[i]=="1|1") {data[i]<-2}
  }
  return(data)
}
fun_conta_ALT_2<-function(data, lista_pop){
  ##################  prendo gli estremi di ogni popolazione per poter poi calcolare i dati mancanti per pop per filtrare dopo
  born_inf_fin<-c()
  born_sup_fin<-c()
  n_pop<-length(lista_pop)
  for (z in 1:n_pop){
    ####ciclo per prendere i limiti di ciascuna pop dove vedere se lo SNP c'? in almeno un individuo
    if (z==1) {
      born_inf=1
      born_sup=length(lista_pop[[z]])
    }
    else {
      vett_lungh_sup<-c()
      for (s in 1:z) {
        vett_lungh_sup<-cbind(vett_lungh_sup,length(lista_pop[[s]]))}
      vett_lungh_inf<-c()
      for (s in 1:(z-1)) {vett_lungh_inf<-cbind(vett_lungh_inf,length(lista_pop[[s]]))}
      
      born_inf=1+sum(vett_lungh_inf)
      born_sup=sum(vett_lungh_sup)
    }
    born_inf_fin<-c(born_inf_fin,born_inf) ## i due vettori con gli indici degli estremi delle popolazioni
    born_sup_fin<-c(born_sup_fin,born_sup)
  }
  ####################################
  risultato<-c()
  #for (i in 1:length(data)) {
  for (j in 1:length(born_inf_fin)) {
    risultato<-c(risultato, sum(data[born_inf_fin[j]:born_sup_fin[j]], na.rm=T))
  }
  risultato<-c(risultato,sum(data, na.rm=T))
  #}
  return(risultato)
}
calcola_sfs2D_pairwise.modified.for.blocks<-function(mat_con_alt,lista_pop){
  n_ind<-2*length(unlist(lista_pop))
  ########## creo matrici di output
  new_list<-list()
  for (i in 1:(length(lista_pop)-1)){
    for (j in 2:length(lista_pop)){
      if (i!=j && i<j) {
        righe=2*length(lista_pop[[i]])+1
        colonne=2*length(lista_pop[[j]])+1
        matrice<-matrix(rep(0,righe*colonne),nrow=righe, ncol=colonne)
        #print(i)
        #print(j)
        new_list<-list.append(new_list,matrice)
      }
    }
  }
  ########
  if(dim(mat_con_alt)[1]!=0){
    #pb <- txtProgressBar(min = 0, max = nrow(mat_con_alt), style = 3)
    for (z in 1:nrow(mat_con_alt)){
      if (mat_con_alt[z,ncol(mat_con_alt)]<(n_ind/2)) { ### il MAF ? l'allele ALT
        indice_list<-1
        for (i in 1:(length(lista_pop)-1)){
          for (j in 2:length(lista_pop)){
            if (i!=j && i<j) {
              new_list[[indice_list]][mat_con_alt[z,i]+1,mat_con_alt[z,j]+1]<-new_list[[indice_list]][mat_con_alt[z,i]+1,mat_con_alt[z,j]+1]+1
              indice_list<-indice_list+1
            }
          }
        }
      } else if (mat_con_alt[z,ncol(mat_con_alt)]>(n_ind/2)) { ### il MAF ? l"allele REF
        indice_list<-1
        for (i in 1:(length(lista_pop)-1)){
          for (j in 2:length(lista_pop)){
            if (i!=j && i<j) {
              n_ind_pop_row<-2*length(lista_pop[[i]])
              n_ind_pop_col<-2*length(lista_pop[[j]])
              new_list[[indice_list]][n_ind_pop_row-mat_con_alt[z,i]+1,n_ind_pop_col-mat_con_alt[z,j]+1]<-new_list[[indice_list]][n_ind_pop_row-mat_con_alt[z,i]+1,n_ind_pop_col-mat_con_alt[z,j]+1]+1
              indice_list<-indice_list+1
            }
          }
        }
      } else {  ### i due alleli hanno freq 0.5 nella pop totale, quindi seguo Laurent che aumenta di 0.5 la entries per entrambi
        indice_list<-1
        for (i in 1:(length(lista_pop)-1)){
          for (j in 2:length(lista_pop)){
            if (i!=j && i<j) {
              n_ind_pop_row<-2*length(lista_pop[[i]])
              n_ind_pop_col<-2*length(lista_pop[[j]])
              new_list[[indice_list]][n_ind_pop_row-mat_con_alt[z,i]+1,n_ind_pop_col-mat_con_alt[z,j]+1]<-new_list[[indice_list]][n_ind_pop_row-mat_con_alt[z,i]+1,n_ind_pop_col-mat_con_alt[z,j]+1]+0.5
              new_list[[indice_list]][mat_con_alt[z,i]+1,mat_con_alt[z,j]+1]<-new_list[[indice_list]][mat_con_alt[z,i]+1,mat_con_alt[z,j]+1]+0.5
              indice_list<-indice_list+1
            }
          }
        }
      }
      #setTxtProgressBar(pb, z)
    }
  }
  return(new_list)
}
make.blocks=function(dati2,mat_con_alt,block_size=10000){
  POS=getPOS(dati2)
  blocks=c(seq(0,max(POS),block_size),max(POS))
  
  mat_con_alt_blocks = list()
  pb <- txtProgressBar(min = 0, max = (length(blocks)-1), style = 3)
  for(i in 1:(length(blocks)-1)){
    mat_con_alt_blocks[[i]] = mat_con_alt[which(POS > blocks[i] & POS <= blocks[i+1]),]
    if(is.null(dim(mat_con_alt_blocks[[i]]))){
      mat_con_alt_blocks[[i]] = t(matrix(mat_con_alt_blocks[[i]]))
    }
    setTxtProgressBar(pb, i)
  }
  return(mat_con_alt_blocks)
}
get.back.full.sfs = function(res,subset_vect=NULL){
  
  if(is.null(subset_vect)){
    subset_vect = c(1:length(res))
  }
  
  res_all = list()
  res_temp = res[[subset_vect[1]]]
  for(j in 1:length(res_temp)){
    res_all[[j]] = res_temp[[j]]
  }
  
  #pb <- txtProgressBar(min = 0, max = (length(blocks)-1), style = 3)
  for(i in 2:length(subset_vect)){
    res_temp = res[[subset_vect[i]]]
    for(j in 1:length(res_temp)){
      res_all[[j]] = res_all[[j]] + res_temp[[j]]
    }
    #setTxtProgressBar(pb, i)
  }
  return(res_all)
}
## Wrapping function
bootSFS2D.blocks=function(dati,lista_pop,block_size=10000,n_boot=100,return_OBS=T){
  
  # This is a part with the function of Stefano
  n_ind<-2*length(unlist(lista_pop)) ### n? of chromosomes, to know who is the MAF
  dati2<-dati[,c("FORMAT", unlist(lista_pop))] ### riordino individui per assegnare pop
  genotype<-extract.gt(dati2, element = "GT", mask = FALSE, as.numeric=F,return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE) ### prendo i genotipi in 1/0
  cat('\n',"Converting genotype file",'\n')
  qqq=pbapply(genotype, 1, fun_geno_mod)
  genotype_num<-matrix(as.numeric(qqq), ncol=ncol(genotype), nrow=nrow(genotype), byrow=T) ### a questo devo applicare funziona conta ALT
  cat('\n',"Computing allele frequencies",'\n')
  mat_con_alt<-t(pbapply(genotype_num, 1, fun_conta_ALT_2, lista_pop)) ### each column is the sim of ALT for each pop. Last colomn is the total ALT, which is needed to know who is the global MAF
  
  
  cat('\n',"Splitting the dataset into blocks of",block_size,"bp",'\n')
  mat_con_alt_blocks = make.blocks(dati2,mat_con_alt,block_size = block_size)
  
  cat('\n',"Computing 2DSFS for each block",'\n')
  res = pbmapply(mat_con_alt_blocks,FUN = calcola_sfs2D_pairwise.modified.for.blocks,MoreArgs=list(lista_pop),SIMPLIFY = F)
  
  ### Just to check if it worked
  #res_all = get.back.full.sfs(res)
  
  #test = calcola_sfs2D_pairwise.mod_boot(mat_con_alt2)
  
  #res_all[[1]] == test[[1]]
  ### Ok we're able to get back the same SFS. 
  
  
  ### Now we bootstrap the sfs 
  cat('\n',"Computing",n_boot,"bootstraps on the 2DSFS",'\n')
  n_boot = n_boot
  boot_sfs = list()
  pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
  for(i in 1:n_boot){
    new_sample = sample(length(res),replace = T)
    new_sample = new_sample[order(new_sample)]
    boot_sfs[[i]] = get.back.full.sfs(res = res,subset_vect = new_sample)
    setTxtProgressBar(pb, i)
  }
  if(return_OBS==T){
    obs_SFS=get.back.full.sfs(res = res,subset_vect = NULL)
    final = list(obs_SFS,boot_sfs); names(final) = c("Observed","Bootstrapped")
  }else{
    final = boot_sfs
  }
  return(final)
} 

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





# ===================================== #
# >>>>>> II. Stairway Plot Tools <<<<<< #
# ===================================== #


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

fluctuations.stairway <- function(stairway,N_DISCRETIZE){
  require(arules)
  
  disc = discretize(stairway$mutation_per_site,breaks = N_DISCRETIZE)
  
  lev = levels(disc)
  
  disc = data.frame(disc)
  
  X = Y = c()
  for(d in 1:length(lev)){
    tempx = stairway$mutation_per_site[which(disc[,1]==lev[d])]
    tempy = stairway$theta_per_site_median[which(disc[,1]==lev[d])]
    X = c(X,mean(tempx))
    Y = c(Y,mean(tempy))
  }
  
  vect_slope = c()
  for(j in 2:length(X)){
    x1 = X[j-1]
    x2 = X[j]
    y1 = Y[j-1]
    y2 = Y[j]
    
    slope = sqrt(((y2-y1)/(x2-x1))^2)
    
    vect_slope <-c(vect_slope,slope)
  }
  return(sum(vect_slope,na.rm = T))
}


# ========================================== #
# >>>> III. FST and individual distance <<<< #
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

genetic.distance=function(gt){
  res = c()
  for(j in 1:length(gt)){
    for(k in j:length(gt)){
      temp = gt[c(j,k)]
      if(length(which(temp=="./." | is.na(temp)))>0){
        res=c(res,NA)
      }else{
        m=matrix(unlist(strsplit(temp,"[/|]")),byrow = T,ncol=2)
        l=c()
        for(i in 1:length(m[1,])){
          if(m[1,1]!=m[1,2]){
            l=c(l,length(which(m[2,] == m[1,i])))
          }else{
            l=length(which(m[2,] == m[1,i]))
          }
        }
        shared_alleles = sum(l)
        d=1-(shared_alleles/2)
        res=c(res,d)
      }
    }
  }
  return(res)
}

BC.distance.individual <- function(vcf,bootstrap=NULL){
  
  GT<-extract.gt(vcf, element = "GT", mask = FALSE,
                 as.numeric=F,return.alleles = FALSE, 
                 IDtoRowNames = TRUE, extract = TRUE, convertNA = T)
  
  BC<-function(GT){
    
    allele.frequency<-function(gt){
      sample_size<-length(gt)*2 - length(which(gt=="./." | gt=="." | gt==".|." | is.na(gt)))*2
      ref<-length(which(unlist(strsplit(gt,"[/|]"))=="0")) / sample_size
      #der<-length(which(unlist(strsplit(gt,"[/|]"))=="1")) / sample_size
      return(ref)
    }
    
    g1<-matrix(GT[,1],ncol = 1)
    g2<-matrix(GT[,2],ncol = 1)
    alleles_frequency1 <- as.numeric(apply(g1,1,allele.frequency))
    alleles_frequency2 <- as.numeric(apply(g2,1,allele.frequency))
    
    AF<-cbind(alleles_frequency1,alleles_frequency2)
    
    AF<-AF[-which(is.na(AF[,1]) | is.na(AF[,2])),]
    
    BC<-(1-2*(sum(apply(AF,1,min)) / (sum(AF[,1]) + sum(AF[,2]))))
    
    return(BC)
  }
  
  resample.alleles<-function(gt){
    if(length(which(is.na(gt))>0)){
      return(c(NA,NA))
    }else{
      new<-sample(size = 2,4)
      g1<-paste0(unlist(strsplit(gt,"[/|]"))[new][1],"/",unlist(strsplit(gt,"[/|]"))[new][1])
      g2<-paste0(unlist(strsplit(gt,"[/|]"))[-new][1],"/",unlist(strsplit(gt,"[/|]"))[-new][1])
      return(c(g1,g2))
    }
    
  }
  
  if(is.null(bootstrap)==FALSE){
    boot_bc<-c()
    for(j in 1:bootstrap){
      temp<-t(apply(gt,1,resample.alleles))
      boot_bc<-c(boot_bc,BC(temp))
    }
  }
  
  K<-ncol(GT)
  res<-matrix(NA,ncol=K,nrow=K)
  colnames(res) = rownames(res) = colnames(GT)
  pb <- txtProgressBar(min = 0, max = K-1, style = 3)
  for(i in 1:(K-1)){
    #cat(">>>",i,"\n")
    for(j in i:K){
      res[i,j] = res[j,i] = BC(GT[,c(i,j)])
    }
    setTxtProgressBar(pb, i)
  }
  return(as.dist(res))
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
# >>>> IV. Filtering and extracting info <<< #
# ========================================== #

# > 1. missing data
missing.data<-function(GT=NULL,vcf=NULL, SAMPLE=T, SNPs=T){
  miss.data<-function(gt){
    return(length(which(gt=="./." | is.na(gt) | gt==".|."))/length(gt))
  }
  
  if(is.null(GT)){
    require(vcfR)
    GT<-extract.gt(vcf, element = "GT", mask = FALSE,
                   as.numeric=F,return.alleles = FALSE, 
                   IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  }
  if(SAMPLE==T){
    missing_data_ind <- apply(GT,FUN = miss.data,MARGIN = 2)
  }
  if(SNPs==T){
    missing_data_site <- apply(GT,FUN = miss.data,MARGIN = 1)
  }
  
  if(SAMPLE==T & SNPs==T){
    final <- list(missing_data_site, missing_data_ind)
    names(final) <- c("missing_data_per_site","missing_data_per_sample")
    return(final)
  }else{
    if(SAMPLE==T){return(missing_data_ind)}
    if(SNPs==T){return(missing_data_site)}
  }
}

# > 2. depth
depth <- function(DP=NULL,vcf=NULL, SAMPLE=T, SNPs=T){
  if(is.null(DP)){
    require(vcfR)
    
    DP<-extract.gt(vcf, element = "DP", mask = FALSE,
                   as.numeric=T,return.alleles = FALSE, 
                   IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  }
  if(SAMPLE==T){
    depth_individual<-apply(DP,MARGIN = 2,FUN = mean)
  }
  if(SNPs==T){
    depth_site<-apply(DP,MARGIN = 1,FUN = mean)
  }
  
  if(SAMPLE==T & SNPs==T){
    final <- list(depth_site, depth_individual)
    names(final) <- c("mean_depth_per_site","mean_depth_per_sample")
    return(final)
  }else{
    if(SAMPLE==T){return(depth_individual)}
    if(SNPs==T){return(depth_site)}
  }
}

# > 3. gc content
gc.content <- function(alleles=NULL, vcf=NULL){
  if(is.null(alleles)){
    require(vcfR)
    
    alleles<-extract.gt(vcf, element = "GT", mask = FALSE,
                        as.numeric=F,return.alleles = TRUE, 
                        IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  }
  gc<-function(alleles){
    a<-unlist(strsplit(alleles,split = "/|"))
    a <- a[-which(a=="|" | a=="" | a==".")]
    gc<-length(which(a=="G" | a=="C"))/length(a)
    at<-length(which(a=="A" | a=="T"))/length(a)
    a<-c(gc,at); names(a)=c("GC","AT")
    return(a)
  }
  
  content<-apply(alleles,MARGIN = 2,FUN = gc)
  return(content)
}

# > 4. proportion of heterozygotes
heterozygosity <- function(GT=NULL,vcf=NULL, SAMPLE=F, SNPs=T){
  if(is.null(GT)){
    require(vcfR)
    require(pbapply)
    
    GT<-extract.gt(vcf, element = "GT", mask = FALSE,
                   as.numeric=F,return.alleles = FALSE, 
                   IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  }
  
  hetero<-function(gt){
    Het<-length(which(gt=="1/0" | gt=="0/1" | gt=="0|1" | gt=="1|0"))
    Total <- length(gt)-length(which(is.na(gt) | gt == "./."))
    return(Het/Total)
  }
  
  
  if(SAMPLE==T){
    hetero_samples<-pbapply(GT,MARGIN = 2,FUN = hetero)
  }
  if(SNPs==T){
    hetero_sites<-pbapply(GT,MARGIN = 1,FUN = hetero)
  }
  
  if(SAMPLE==T & SNPs==T){
    final <- list(hetero_sites, hetero_samples)
    names(final) <- c("heterozygosity_per_site","heterozygosity_per_sample")
    return(final)
  }else{
    if(SAMPLE==T){return(hetero_samples)}
    if(SNPs==T){return(hetero_sites)}
  }
}

# > 5. proportion of fixed alleles
fixed.alleles <- function(GT=NULL,vcf=NULL, SAMPLE=F, SNPs=T){
  if(is.null(GT)){
    require(vcfR)
    
    GT<-extract.gt(vcf, element = "GT", mask = FALSE,
                   as.numeric=F,return.alleles = FALSE, 
                   IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  }
  fixed.ref<-function(gt){
    return(length(which(gt=="0/0" | gt=="0|0"))/length(gt))
  }
  fixed.der<-function(gt){
    return(length(which(gt=="1/1" | gt=="1|1"))/length(gt))
  }
  if(SAMPLE==T){
    prop_ref_ind<- apply(GT,FUN = fixed.ref,MARGIN = 2)
    prop_der_ind<- apply(GT,FUN = fixed.der,MARGIN = 2)
  }
  if(SNPs==T){
    prop_ref_site <- apply(GT,FUN = fixed.ref,MARGIN = 1)
    prop_der_site <- apply(GT,FUN = fixed.der,MARGIN = 1)
  }
  
  if(SAMPLE==T & SNPs==T){
    final <- list(prop_ref_ind, prop_der_ind, prop_ref_site, prop_der_site)
    names(final) <- c("ref_sample","der_sample","ref_site","der_site")
    
  }else{
    if(SAMPLE==T){
      final <- list(prop_ref_ind, prop_der_ind)
      names(final) <- c("ref_sample","der_sample")
      
    }
    if(SNPs==T){
      final <- list(prop_ref_site, prop_der_site)
      names(final) <- c("ref_site","der_site")
      
    }
  }
  return(final)
}

# > 6. allele frequency
minor.allele.frequency<-function(GT=NULL,vcf=NULL){
  require(pbapply)
  if(is.null(GT)){
    require(vcfR)
    GT<-extract.gt(vcf, element = "GT", mask = FALSE,
                   as.numeric=F,return.alleles = FALSE, 
                   IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  }
  
  
  Minor.Allele.Frequency<-function(gt){
    alleles <- unlist(strsplit(gt,split = "[/|]"))
    A <- length(which(alleles=="0"))/length(alleles)
    B <- length(which(alleles=="1"))/length(alleles)
    return(min(A,B))
  }
  
  maf<-pbapply(GT,MARGIN = 1,FUN = Minor.Allele.Frequency)
  return(maf)
}
multi.minor.allele.frequency<-function(GT=NULL,vcf=NULL){
  require(pbapply)
  if(is.null(GT)){
    require(vcfR)
    GT<-extract.gt(vcf, element = "GT", mask = FALSE,
                   as.numeric=F,return.alleles = FALSE, 
                   IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  }
  
  #gt=GT[19,]
  
  Multi.Minor.Allele.Frequency<-function(gt){
    alleles <- unlist(strsplit(gt,split = "[/|]"))
    REF <- length(which(alleles=="0"))/length(alleles)
    
    all = length(levels(as.factor(alleles)))
    ALT = c()
    for(i in 1:(all-1)){
      ALT <- c(ALT,length(which(alleles==as.character(i)))/length(alleles))
    }
    
    return(min(REF,ALT))
  }
  
  maf<-pbapply(GT,MARGIN = 1,FUN = Multi.Minor.Allele.Frequency)
  return(maf)
}
allele.frequecy.ref <- function(GT=NULL,vcf=NULL){
  require(pbapply)
  
  if(is.null(GT)){
    require(vcfR)
    GT<-extract.gt(vcf, element = "GT", mask = FALSE,
                   as.numeric=F,return.alleles = T, 
                   IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  }
  
  allREF = getFIX(vcf)[,"REF"] 
  
  GT = cbind(GT,allREF)
  
  AFR=function(gt){
    ref = gt["allREF"]
    gt = gt[1:(length(gt)-1)]
    alleles <- unlist(strsplit(gt,split = "[/|]"))
    return(length(which(alleles==ref)) / length(alleles))
  }
  
  REF = pbapply(GT,MARGIN = 1,FUN = AFR)
  
  #REF <- c()
  #pb <- txtProgressBar(min = 0, max = nrow(GT), style = 3)
  #for(i in 1:nrow(GT)){
  #  alleles <- unlist(strsplit(GT[i,],split = "[/|]"))
  #  REF = c(REF,length(which(alleles==allREF[i])) / length(alleles))
  #  setTxtProgressBar(pb, i)
  #}
  n = length(unlist(strsplit(rownames(GT)[1],"_")))
  ssss =as.numeric(matrix(unlist(strsplit(rownames(GT),"_")),ncol = n,byrow = T)[,3])
  return(cbind(ssss,REF))
}

# > 7. pick random SNP per locus
random.snp<-function(vcf){
  require(vcfR)
  GT<-extract.gt(vcf, element = "GT", mask = FALSE,
                 as.numeric=F,return.alleles = FALSE, 
                 IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)
  
  loci<-getCHROM(vcf)
  loci_unique<-loci[-which(duplicated(loci))]
  
  to_remove<-c()
  for(i in 1:length(loci_unique)){
    if(length(which(loci==loci_unique[i]))>1){
      snps<-which(loci==loci_unique[i])
      temp<-sample(x = snps,size = (length(snps) - 1) ,replace = F)
      to_remove<-c(to_remove,temp)
    }
    p<-round((i/length(loci_unique))*100)
    if(p!=round(((i-1)/length(loci_unique))*100)){
      cat("=")
    }
  }
  
  length(to_remove)
  
  vcf_clean<-vcf[-to_remove,]
  return(vcf_clean)
}

# > 8. Make a sliding window coordinate dataset
sliding.window.df = function(df, jump,slide){
  require(pbapply)
  AF=df$AF
  positions<-df$POS
  nw<-data.frame(low_bound=seq(1, max(positions) + jump, jump),
                 upper_bound=seq(slide, max(positions)+jump+slide, jump))
  nw.list <- split(nw, seq(nrow(nw)))
  sliding.avg<-function(wind,positions){
    temp_res<-matrix(NA,nrow = 1)
    temp<-which(positions >= as.numeric(wind[1]) & positions < as.numeric(wind[2]))
    if(length(temp)>0){
      q_temp<-AF[temp]
      res=mean(q_temp)
    }else{res=NA}
    return(res)
  }
  res = pbmapply(sliding.avg,nw.list,MoreArgs = list(positions))
  res = data.frame(POS=c((nw[,2]+nw[,1])/2),AVG_AF=res)
  return(res)
}

# > 9. Finding discriminating SNPs between two groups of individuals
len.fix.sites <- function(data){
  
  fix <- function(af1, af2, positions,all=0){
    res <- c()
    pb <- txtProgressBar(min = 0, max = length(af1), style = 3)
    for(i in 1:length(af1)){
      if(i != 1 & i != length(af1)){
        if(af1[i] == all){
          k=0
          other_keep = posikeep = c()
          while(af1[i+k] == all){
            other_keep<-c(other_keep,af2[i+k])
            posikeep = c(posikeep,positions[i+k])
            k = k + 1
          }
          STOP=(positions[i+k] - 1)
          START=(positions[i-1] + 1)
          len =  STOP - START 
          SNPs2 = length(which(other_keep!=0 & other_keep!=1))
          temp = data.frame(START=START, STOP=STOP, length=len,
                            N_Fixed=k, N_SNPs_alt_cl=SNPs2)
          res = rbind(res,temp)
        }
      }
      setTxtProgressBar(pb, i)
    }
    res = res[-which(duplicated(res$STOP)),]
    return(res)
  }
  
  cat("\n","1. computing fixed alleles for ref in cluster 1, SNPs in cluster 2","\n")
  fix1 = fix(af1 = data[,2],af2 = data[,3],positions = data[,1],all = 0)
  #cat("\n","2. computing fixed alleles for ref in cluster 2, SNPs in cluster 1","\n")
  #fix2 = fix(af1 = data[,3],af2 = data[,2],positions = data[,1],all = 0)
  cat("\n","3. computing fixed alleles for ALT in cluster 1, SNPs in cluster 2","\n")
  fix1b = fix(af1 = data[,2],af2 = data[,3],positions = data[,1],all = 1)
  #cat("\n","4. computing fixed alleles for ALT in cluster 2, SNPs in cluster 1","\n")
  #fix2b = fix(af1 = data[,3],af2 = data[,2],positions = data[,1],all = 1)
  
  return(list(fix1,fix1b))
}



