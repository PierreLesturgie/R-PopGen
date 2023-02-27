#################################################################
#---------------------------------------------------------------#
#------------- FUNCTIONS TO EXTRACT INFO FROM VCF --------------#
#----------------------------------------------------------------#
#################################################################

require(vcfR)

# > 1. missing data
missing.data<-function(GT=NULL,vcf=NULL, SAMPLE=T, SNPs=T){
  miss.data<-function(gt){
    return(length(which(gt=="./." | is.na(gt)))/length(gt))
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


# > 2. mean depth
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

# > 6. minor allele frequency
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

select.indep.snps<-function(bed_file,chr=c(paste0("SUPER_",1:10)),sep_len=1000){
  if(is.character(bed_file)==F){
    stop("Provide pathway and name of the bed file")
  }
  
  bed<-read.table(bed_file)[,1:2]
  bed[,2] <- bed[,2] + 1 #### I just transform the bed (start 0) in a "position" file useful for after
  
  #==> FIRST WE SPLIT
  bed_list<-list()
  for(i in 1:length(chr)){
    bed_list[[i]]<-bed[which(bed[,1]==chr[i]),]
  }
  
  if(is.null(regions_to_remove)==F){
    names(bed_list) = names(regions_to_remove) = chr
  }else{names(bed_list) = chr}
  
  
  for(i in 1:length(chr)){
    if(is.null(regions_to_remove)==F){
      temp_region <- regions_to_remove[[i]]
      temp_chr <- bed_list[[i]]
      
      idx_regions <- c() ### removing regions
      for(j in 1:nrow(temp_region)){
        idx_regions <- c(idx_regions,which(temp_chr[,2]>=temp_region[j,1] & temp_chr[,2]<=temp_region[j,2]))
      }
      if(length(idx_regions)>0){
        temp_chr <- temp_chr[-idx_regions,]
      }
      cat(">>> Removed",length(idx_regions), "out of", nrow(temp_chr),"SNPs for chromosome", names(bed_list)[i],"\n")
    }else{temp_chr <- bed_list[[i]]}
    
    
    ### choosing first SNP position to consider
    temp_first <- which(temp_chr[,2] >= temp_chr[1,2] & temp_chr[,2] < (temp_chr[1,2]+sep_len))
    keep<-sample(temp_first,size = 1)
    pos<-temp_chr[keep,2]
    
    cat(">>> Starting at position",pos, "for chromosome", names(bed_list)[i],"\n")
    idx_pos<-c(1:(keep-1))
    for(l in (keep+1):nrow(temp_chr)){
      if((temp_chr[l,2]-pos) < sep_len){
        idx_pos<-c(idx_pos,l)
      }else{pos<-temp_chr[l,2]}
      # cat(chr[i],l,"/",nrow(temp_chr),"\n")
    }
    temp_chr <- temp_chr[-idx_pos,]
    
    cat(">>> Removed",length(idx_pos), "out of", l,"SNPs for chromosome", names(bed_list)[i],"\n")
    
    write.table(temp_chr,paste0("filtered_",chr[i],".txt"),
                sep = "\t",quote = F,row.names = F,col.names = F)
    
  }
}