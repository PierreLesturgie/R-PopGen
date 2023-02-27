#################################################################
#---------------------------------------------------------------#
#------------------------ Hudson's FST -------------------------#
#---------------------------------------------------------------#
#################################################################



# 1. translate genotypes (as in extracted by VCFR) to alleles
gt2alleles = function(gt){
  alleles = unlist(strsplit(gt,"[/|]"))
  alleles[which(alleles==".")] = NA
  while(length(alleles) != length(gt)*2){
    alleles <- c(alleles,NA)
  }
  return(as.numeric(alleles))
}

# 2. compute mean pairwise difference from allele number and allele counts
mpd = function(an,ac){
  n_pairs = an * (an - 1) / 2
  n_same = rowSums(ac * (ac - 1) / 2)
  
  n_diff = n_pairs - n_same
  
  return(n_diff / n_pairs)
}

# 3. compute mean pairwise difference from allele number and allele counts between populations
mpd.between = function(an1,an2,ac1,ac2){
  n_pairs = an1 * an2 
  n_same = rowSums(ac1 * ac2)
  
  n_diff = n_pairs - n_same
  
  return(n_diff / n_pairs)
}

# 4. Computes a pairwise FST from allele counts of the two populations
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

# 4. Computes a pairwise FST from allele counts of the two populations in a sliding window fashion
sliding.window.fst = function(windows,positions,an1,an2,ac1,ac2){
  temp<-which(positions >= as.numeric(windows[1]) & positions < as.numeric(windows[2]))
  if(length(temp) > 0){
    pairwsie_fst = p.FST(an1 = an1[temp], an2 = an2[temp], 
                         ac1 = ac1[temp,], ac2 = ac2[temp,])
    return(pairwsie_fst[[1]])
  }else{return(NA)}
}

# 5. calculate number of alleles and filters 
a.counts = function(site,mx){
  c <- c()
  for(i in 1:(mx+1)){
    c <- c(c,length(which(site == (i-1))))
  }
  return(c)
}

# 6. Wrapping function - only this one should be used  
# >>> Only required arguments are a vcf and a poplist
# >>> The function implements resampling - but takes time. 
fst.hudson <- function(vcf, pop_list,resampling=NULL,sliding_window=FALSE,slide=NULL,jump=NULL,write=FALSE){
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

# 7. creates a distance matrix from geo coordinates (lon, lat)
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

# 8. Performs a Mantel test to test isolation by distance
Isolation.By.Distance<-function(dgeo,dfst){
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

# 9. Tools to plot FST: outputs heatmap and network
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