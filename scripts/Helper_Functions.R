## Function to cluster the data
.cluster_dat <- function(x, metric, linkage, rcpp) {
   if (metric == "cor") {
    dmat <- 1 - WGCNA::cor(t(x))
    hc_dat <- hclust(as.dist(dmat), method=linkage)
  } else if (rcpp) {
    hc_dat <- Rclusterpp::Rclusterpp.hclust(x, method=linkage, 
                                            distance=metric)
  } else {
    hc_dat <- hclust(dist(x, method=metric), method=linkage)
  }
  hc_dat
}


## perform hierarchical clustering on the original data and 
## compute the corresponding cluster indices for each merge
.initcluster <- function(x, n, p, metric, linkage, 
                        n_ci, ci, rcpp) {
  
  ## obtain clustering solution
  hc_dat <- .cluster_dat(x, metric, linkage, rcpp)
  
  ## list array of cluster indices at each of the n-1 nodes
  idx_hc <- .idx_hc(hc_dat, n)
  
  ## matrix containing cluster indices
  ci_dat <- matrix(-1, nrow=n-1, ncol=n_ci)
  
  ## calculate cluster index(ices) for merge k
  for (i_ci in 1:n_ci) {
    if (ci[i_ci] == "2CI") {
      for (k in 1:(n-1)) {
        ci_dat[k, i_ci] <- .calc2CI(x[idx_hc[[k, 1]], , drop=FALSE],
                                    x[idx_hc[[k, 2]], , drop=FALSE])
      }
    } else if (ci[i_ci] == "linkage") {
      ## flip index, hclust and shc use revered ordering
      ci_dat[, i_ci] <- rev(hc_dat$height)
    }
  }
  
  list(hc_dat = hc_dat, 
       idx_hc = idx_hc,
       ci_dat = ci_dat)
}


## determine obs indices at each node of the dendrogram
.idx_hc <- function(hc, n) {
  ## list array of cluster indices at each of the n-1 merges
  idx_hc <- array(list(), c(2*n-1, 2))
  idx_hc[1:n, 1] <- as.list(n:1)
  idx_hc[(n+1):(2*n-1), ] <- hc$merge + n + (hc$merge<0)
  
  ## complete idx_hc
  for (k in 1:(n-1)) {
    idx_hc[[n+k, 1]] <- unlist(idx_hc[idx_hc[[n+k, 1]], ])
    idx_hc[[n+k, 2]] <- unlist(idx_hc[idx_hc[[n+k, 2]], ])
  }
  
  ## flip index, hclust and shc use revered ordering
  idx_hc[(2*n-1):(n+1), ]
}


## calculate 2-means cluster index (n x p matrices)
.calc2CI <- function(x1, x2) {
  if (is.matrix(x1) && is.matrix(x2) && ncol(x1) == ncol(x2)) {
    (.sumsq(x1) + .sumsq(x2)) / .sumsq(rbind(x1, x2))
  } else {
    stop(paste("x1, x2 must be matrices with same ncols",
               "for 2CI calculation"))
  }      
}

## calculate sum of squares
.sumsq <- function(x) { norm(sweep(x, 2, colMeans(x), "-"), "F")^2 }


## Function to compute differential expression
.Defexp <- function(cluster_1, cluster_2, X, k){
  
  # Compute the means of the samples in each cluster
  c1.mean = apply(cluster_1, 1, mean)
  c2.mean = apply(cluster_2, 1, mean)

  # Just get the maximum of all the means
  limit = max(c1.mean, c2.mean)
  
  # Compute fold-change (biological significance)
  # Difference between the means of the conditions
  fold = c1.mean - c2.mean
  
  
  # Compute statistical significance (using t-test)
  pvalue = NULL # Empty list for the p-values
  tstat = NULL # Empty list of the t test statistics
  
  for(i in 1 : nrow(X)) { # For each gene : 
    x = cluster_1[i,] # WT of gene number i

    y = cluster_2[i,] # KO of gene number i
    
    # Compute t-test between the two conditions
    t = t.test(x, y)
    
    # Put the current p-value in the pvalues list
    pvalue[i] = t$p.value
    # Put the current t-statistic in the tstats list
    tstat[i] = t$statistic
  }
  
  # Screen for the genes that satisfy the filtering criteria
  fold_cutoff = 2

  pvalue_cutoff = 0.01
  # Fold-change filter for "biological" significance
  filter_by_fold = abs(fold) >= fold_cutoff
  
  # P-value filter for "statistical" significance
  filter_by_pvalue = pvalue <= pvalue_cutoff
  X[filter_by_pvalue, ]
  
  # Combined filter (both biological and statistical)
  filter_combined = filter_by_pvalue
  
  filtered = X[filter_combined,]
  
  list("K" = k, "tstat" = tstat, "pvals" = pvalue, "filtered" = filtered)
}


#function to create a seurat object
.CreateSeuratObj <- function(path, Pname, minCells, minGenes, matrix = FALSE, sep = " "){
  #if we are using the 10X matrix
  if(matrix == TRUE){
    data <-  Read10X(data.dir = getwd())
    obj <- CreateSeuratObject(raw.data = data, project = Pname)
    #min.cells = minCells, min.genes = minGenes,
    return(obj)
  }
#read in the data
data <- fread(path) %>% column_to_rownames(var = "V1")
#create the seurat object
obj <- CreateSeuratObject(raw.data = data, project = Pname, min.cells = minCells, min.genes = minGenes)

return(obj)
}

#Create a function to preprocess the data
.PreProcess <- function(obj, low_thres, high_thres){
  #get the mitohondrial genes
  mito.genes <- grep(pattern = "^mt-", x = rownames(x = obj@data), value = TRUE)
  percent.mito <- Matrix::colSums(obj@raw.data[mito.genes, ])/Matrix::colSums(obj@raw.data)
  obj <- AddMetaData(object = obj, metadata = percent.mito, col.name = "percent.mito")
  
  #do some filtering
  obj <- FilterCells(object = obj, subset.names = c("nGene", "percent.mito", "nUMI"), 
                     low.thresholds = low_thres, high.thresholds = high_thres)

  #normalize and scale the data
  obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  #compute highly variable genes
  obj <- FindVariableGenes(object = obj, do.cpp = T,
                           x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 1, do.plot = F)
  #mean.function = ExpMean, dispersion.function = LogVMR,
  
  #scale data and remove unwanted sources of variation
  #obj <- ScaleData(object = obj, vars.to.regress = c("percent.mito", "nUMI"))
  
  return(obj)
}

#create a function to extract the data of the top n highly variable genes
.GetGeneData <- function(obj, top_N, scaled = T){
  #extract the data
  if(scaled == T){
    obj_sub <- as.matrix(obj@scale.data[top_N, ])
  }
  if(scaled == F){
    obj_sub <- as.matrix(obj@data[top_N, ])
  }
  
  return(obj_sub)
}


#create a function to sort the variable genes by dispersion
.SortGenes <- function(obj){
  
  #extract the variable genes from the seurat object
  var_genes <- obj@hvg.info[obj@var.genes, ]
  
  #add a column for the genes
  var_genes <- setDT(var_genes, keep.rownames = TRUE)[]
  colnames(var_genes) <- c("gene", "mean", "dispersion", "scaled.dispersion")
  
  #sort by scaled dispersion in descending order
  var_genes <- var_genes[order(-var_genes$scaled.dispersion), ]
  
  #return the results
  return(var_genes)
}


#create a function to run the pipeline
.RunPipeline <- function(Ngenes, obj, sig_level = 0.05, De_plot = TRUE){
  #sort the genes
  var_genes_sorted <- .SortGenes(obj)
  
  #select the top N genes
  top_N <- var_genes_sorted$gene[1:Ngenes]
  
  #select the top N genes
  obj_sub <- .GetGeneData(obj, top_N)
  
  #Run significant clustering
  obj_shc <- shc(t(obj_sub), matmet = .fast_dist, linkage="ward", alpha=0.05, rcpp = T, n_min = 50)
  p <- plot(obj_shc, hang=.1, main = 'obj')
  
  #compute the number of significant and non significant nodes
  obj_sub <- t(obj_sub)
  sig = obj_shc$nd_type == "sig"
  sig_ind <- obj_shc$idx_hc[sig, ]
  
  #get the indices of the significant nodes
  sig_arr <- which(sig, arr.ind = TRUE, useNames = TRUE)
  not_sig <- obj_shc$nd_type == "not_sig"
  not_sig_ind <- obj_shc$idx_hc[not_sig, ]
  
  #perform differential expression for both significant and non significant nodes
  sig_res <- NULL
  not_sig_res <- NULL
  genes <- top_N
  
  print("Now running differential expression for significant nodes")
  #perform DE for significant nodes
  DE_arr = NULL
  DE_genes <- NULL
  if(length(sig_ind) > 0){
    count <- 1
    x = length(sig_ind)/2
    for (i in 1:x){
      c1 <- t(obj_sub[sig_ind[[i]], ])
      c2 <- t(obj_sub[sig_ind[[i+x]], ])
      Xnew <- cbind(c1, c2)
      
      #run DE
      res <- DiffTTest(obj, cells.1 = colnames(c1), cells.2 = colnames(c2), genes.use = obj@var.genes)
      pvals = p.adjust(res$p_val, method = "bonferroni", n = x*length(obj@var.genes))
      res$p_val_adj <- pvals
      #compute the % percentage of DE genes.
      num_sig <- length(res$p_val_adj[res$p_val_adj <= sig_level])
      #extract the significant genes
      temp_genes <- res[res[, "p_val_adj"] <= sig_level,]
      #sort the significant genes
      temp_genes <- temp_genes[order(temp_genes$p_val_adj),]
      temp_genes$genes <- rownames(temp_genes)
      DE_genes[[i]] <- temp_genes
      
      
      DE_arr[[i]] <- num_sig
      print(num_sig)
      sig_res[[i]] <- res
      
      
    }

  }
  
  else{
    return(NULL)
  }
  #update the plot values
  for(i in sig_arr){ 
    obj_shc$p_norm[i] <- DE_arr[i]/10000
    }
  
  if(De_plot == TRUE){
    p_new <- plot(obj_shc, hang=.1, main = 'obj')
    return(list("sig" = sig_res, "plot" = p, "Sig_DE" = DE_arr, "sigclust" = obj_shc, "DE_plot" = p_new, "DE_genes" = DE_genes))
  }
  
  return(list("sig" = sig_res, "plot" = p, "Sig_DE" = DE_arr, "sigclust" = obj_shc, "DE_genes" = DE_genes))
}


#create a function to return the union of all the significant DE genes
ComputeSet <- function(results, type = "union"){
  # a list to store the union values
  set_vals <- NULL
  #loop over all the gene set tested
  for(i in 1:length(results)){
    temp <- results[[i]]$DE_genes
    #check if more than one node was significant for this run.
    if(length(temp) > 1){
      temp_vals <- as.data.frame(temp[[1]]$genes)
      colnames(temp_vals) <- c("genes")
      #loop over all significant nodes
      for(j in 2:length(temp)){
        curr <- as.data.frame(temp[[j]]$genes)
        colnames(curr) <- c("genes")
        if(type == "union"){
          temp_vals <- union(temp_vals, curr)
        }
        if(type == "intersect"){
          temp_vals <- intersect(temp_vals, curr)
        }
      }
      #reorder the columns and update the union list
      #temp_unions <- temp_unions[, c(3,1,2)]
      #temp_unions <- temp_unions[order(temp_unions$p_val_adj), ]
      temp_vals <- as.data.frame(temp_vals)
      
      #remove any NA values
      temp_vals <- temp_vals[!grepl("NA",temp_vals$genes),]
      #temp_unions <- unique( temp_unions[ , c(1) ] )
      if(length(temp_vals) > 0){
        set_vals[[i]] <- temp_vals
      }
      
    }
    else{
      #reorder the columns
      temp <- temp[[1]]$genes
      #temp <- temp[, c(3,1,2)]
      #sort by adjusted pvalues
      #temp <- temp[order(temp$p_val_adj), ]
      temp <- as.data.frame(temp)
      temp <- temp[!grepl("NA",temp$genes),]
      #temp <- unique( temp[ , c(1) ] )
      #update the union list
      set_vals[[i]] <- temp
    }
  }
  #return the union list
  return(set_vals)
}

.GenerateConsensus <- function(Ngenes, obj){
  #sort the genes
  var_genes_sorted <- .SortGenes(obj)
  
  #select the top N genes
  top_N <- var_genes_sorted$gene[1:Ngenes]
  
  #select the top N genes
  obj_sub <- .GetGeneData(obj, top_N)
  
  #compute the distance matrix
  dist_mat <- parallelDist(t(obj_sub), method = "euclidean")
  
  #run the clustering algorithm
  dend <- hclust(dist_mat, method = "ward.D2")
  
  #return the results
  return(dend)
  
}

#write a function to do bootstrap clustering
.BootClust <- function(n_iter = 100, n_genes = 100, dat){
  clust_list <- NULL
  for(i in 1:n_iter){
    
    #sample the data
    permuted <- dat[sample(nrow(dat), n_genes, replace=F),]
    
    #compute the distance matrix
    dist_mat <- parallelDist(t(permuted), method = "euclidean")
    
    #run the clustering algorithm
    dend <- hclust(dist_mat, method = "ward.D2")
    clust_list[[i]] <- dend
  }
  return(clust_list)
}


.Boot_Con <-
  function(data, obj, distance=c("binary","euclidean","maximum","manhattan","canberra",
                           "minkowski","gower","chisq", "correlation"),method=c("complete","ward.D2","single","average","mcquitty","median",
                                                                 "centroid"),nboot=500,duplicate=TRUE,cex.text=1,col.text="red", ...)
  {
    t0<-Sys.time()
    distance <- match.arg(distance)
    method <- match.arg(method)
    if(distance=="gower") {
      if (requireNamespace("cluster", quietly = TRUE)) {
        distancia<-cluster::daisy(data,metric=distance)
      }
      else{
        return("Please install cluster package for calculate gower distance")
      }
    }
    if(distance=="chisq") {
      n<-sum(data)
      nr<-nrow(data);nc<-ncol(data)
      ro<-apply(data,1,function(x)sum(x,na.rm=TRUE))
      co<-apply(data,2,function(x)sum(x,na.rm=TRUE))/n
      B<-data/ro
      A<-matrix(0,nr,nr)
      colnames(A)<-rownames(A)<-rownames(data)
      for(i in 1:nr){
        for(j in 1:nr){d<-0
        for(k in 1:nc) d<-d+(B[i,k]-B[j,k])^2/co[k]
        A[i,j]<-sqrt(d)
        }}
      distancia<-as.dist(A)
    }
    if(distance!="chisq" & distance!="gower" & distance != "correlation"){
      distancia<-parallelDist(data,method=distance)
    }
    
    if(distance == "correlation"){
      distancia <- as.dist(1 - WGCNA::cor1(t(data)))
    }
    
    nc<-ncol(data)
    nr<-nrow(data)
    dend<-fastcluster::hclust(distancia,method=method)
    h1<-cutree(dend,h=0)
    h2<-data.frame(h1)
    h3<-unique(h2)
    dup<-duplicate
    duplicate<-NULL
    # To study duplicate
    if(dup){
      if(nrow(h3) < length(h1)){
        nr<-nrow(h3)
        data<-data.frame(d=rownames(data),data)
        h3<-data.frame(d=rownames(h3),h3)
        duplicate<- merge(h3,data,by="d",all=TRUE)
        duplicate<-duplicate[is.na(duplicate[,2]),]
        dup0<-duplicate[,-2]
        duplicate<-as.character(duplicate$d)
        data<-merge(h3,data,by="d")
        dup1 <-data[,-2]
        dup0<-cbind(dup0,unique="")
        dup0[,1]<-as.character(dup0[,1])
        dup1[,1]<-as.character(dup1[,1])
        ncdup<-ncol(dup1)
        dup0[,ncdup+1]<-as.character(dup0[,ncdup+1])
        ndup0<-nrow(dup0)
        ndup1<-nrow(dup1)
        ncdup<-ncol(dup1)
        for ( i in 1:ndup0) {
          for ( j in 1:ndup1) {
            if(sum(dup0[i,2:ncdup]==dup1[j,-1],na.rm=TRUE) == ncdup-1){
              dup0[i,ncdup+1]<-dup1[j,1]
              break
            }
          }
        }
        if (sum(dup0[,ncdup+1]=="")>0) {
          add1<-dup0[dup0[,ncdup+1]=="",]
          add1<-data.frame(d=add1[,1],h1=0,add1[,2:ncdup])
          data<-rbind(data,add1)
        }
        rownames(data)<-data[,1]
        data<-data[,c(-1,-2)]
        nc<-ncol(data)
        if(distance=="gower") distancia<-daisy(data,metric=distance)
        if(distance=="chisq") {
          n<-sum(data)
          nr<-nrow(data);nc<-ncol(data)
          ro<-apply(data,1,function(x)sum(x,na.rm=TRUE))
          co<-apply(data,2,function(x)sum(x,na.rm=TRUE))/n
          B<-data/ro
          A<-matrix(0,nr,nr)
          colnames(A)<-rownames(A)<-rownames(data)
          for(i in 1:nr){
            for(j in 1:nr){d<-0
            for(k in 1:nc) d<-d+(B[i,k]-B[j,k])^2/co[k]
            A[i,j]<-sqrt(d)
            }}
          distancia<-as.dist(A)
        }
        if(distance!="chisq" & distance!="gower")distancia<-parallelDist(data,method=distance)
        dend<-fastcluster::hclust(distancia,method=method)
        duplicate<-dup0[dup0[,ncdup+1]!="",][,c(1,ncdup+1)]
        names(duplicate)[1]<-"duplicate"
      }}
    nr<-nrow(data)
    if(!is.null(duplicate)) {
      cat("\nDuplicates:",nrow(duplicate))
      cat("\nNew data  :", nr,"Records\n")
    }
    
    dinicio<-dend$merge
    d0<-hgroups(dinicio)
    clases<- data.frame(d=d0,height=dend$height,sq=1:length(d0))
    #rownames(clases)<-d0
    
    b<- nboot
    #b <- seq(50,500, by = 50)
    d<-NULL
    orig_data <- data
    #########################
    for (k in 1:b) {
      #muestra<-sample(1:100,replace=F)
      #sort the genes
      #var_genes_sorted <- .SortGenes(obj)
      
      #select the top N genes
      #top_N <- var_genes_sorted$gene[1:k]
      
      #select the top N genes
      #obj_sub <- .GetGeneData(obj, top_N)
      dat <- orig_data[ ,sample(ncol(orig_data), ceiling(ncol(orig_data)/2), replace=F)]
      #boot1<-data[,muestra]
      boot1 <- t(dat)
      
      data <- t(dat)
      if(distance=="gower") distancia<-daisy(boot1,metric=distance)
      if(distance=="chisq") {
        n<-sum(boot1)
        nr<-nrow(boot1);nc<-ncol(boot1)
        ro<-apply(boot1,1,function(x)sum(x,na.rm=TRUE))
        co<-apply(boot1,2,function(x)sum(x,na.rm=TRUE))/n
        B<-data/ro
        A<-matrix(0,nr,nr)
        colnames(A)<-rownames(A)<-rownames(data)
        for(i in 1:nr){
          for(j in 1:nr){d<-0
          for(k in 1:nc) d<-d+(B[i,k]-B[j,k])^2/co[k]
          A[i,j]<-sqrt(d)
          }}
        distancia<-as.dist(A)
      }
      if(distance!="chisq" & distance!="gower" & distance!="correlation"){
        distancia<-parallelDist(boot1,method=distance)
      }
      if(distance == "correlation"){
        distancia <- as.dist(1 - WGCNA::cor1(t(boot1)))
      }
      d1<-fastcluster::hclust(distancia,method=method)$merge
      d1<-hgroups(d1)
      d<-rbind(d,d1)
    }
    td<-table(d)
    td<-data.frame(td)
    junto<-merge(clases,td,by="d")
    junto[,4]<-junto[,4]*100/b
    junto<-merge(clases,junto,by="d",all=TRUE)
    junto[is.na(junto)]<-0
    junto<-junto[order(junto[,3]),]
    ############
    tiempo<- Sys.time()-t0
    unidad <- attributes(tiempo)$units
    cat("\nConsensus fastcluster::hclust\n" )
    cat("\nMethod distance:",distance)
    cat("\nMethod cluster :",method)
    cat("\nrows and cols  :",nr,nc)
    cat("\nn-bootstrap    :",nboot)
    cat("\nRun time       :",tiempo,unidad,"\n\n")
    cc<-cbind(dend$merge,height=dend$height,percentage=round(junto[,6],1))
    co<-dend$order
    n1<-nrow(cc)
    n2<-length(co)
    p<-rep(0,n1)
    for(i in 1:n1) {
      if ((cc[i,1] < 0) & (cc[i,2] < 0) ) {
        for(k in 1:n2) {
          if(co[k]==-cc[i,1]) k1=k
          if(co[k]==-cc[i,2]) k2=k
        }
        p[i]<- (k1+k2)/2
      }
      if ((cc[i,1]) < 0 & (cc[i,2] > 0)) {
        for(k in 1:n2) {
          if(co[k]==-cc[i,1]) k1=k
        }
        p[i]<-(k1+p[cc[i,2]])/2
      }
      if ((cc[i,1] > 0) & (cc[i,2] < 0)) {
        for(k in 1:n2) {
          if(co[k]==-cc[i,2]) k1=k
        }
        p[i]<-(k1+p[cc[i,1]])/2
      }
      if ((cc[i,1] > 0) & (cc[i,2] > 0)) {
        p[i]<-(p[cc[i,1]]+p[cc[i,2]])/2
      }
    }
    table.dend <- data.frame(dend$merge,xaxis=p,cc[,3:4],groups=d0)
    #plot(dend,...)
    #delta<-0.05*(max(cc[,4]))
    #text(p,cc[,3],ceiling(cc[,4]),cex=cex.text,col=col.text)
    return(list(table.dend=table.dend,dendrogram=dend,duplicates=duplicate) )
  }


## Create a function that takes in the data and returns a list of all the tests

.RunTests <- function(data, iter = 10){
  #Dist dip and silverman tests 
  dist_data <- parallelDist(t(data))
  #DistSilvr <- silverman.test(dist_data,1,adjust=TRUE)@p_value
  DistDip <- dip.test(dist_data)$p.value
  
  #Classic Tests
  #classicSilver <- silverman.test(t(data),1,adjust=TRUE)@p_value
  #classicDip <- dip.test(t(data))$p.value
  
  #PCA tests
  #pcaSilver <-silverman.test(as.matrix(prcomp(t(data))$x[,1]),1,adjust=TRUE)@p_value
  pcaDip <- dip.test(as.matrix(prcomp(t(data))$x[,1]))$p.value
  
  return(list("DipDist" = DistDip, "PCA Dip" = pcaDip))
  #return(list("DipDist" = DistDip, "SilvDist" = DistSilvr, "Dip" = classicDip, "Silv" = classicSilver, "PCA Dip" = pcaDip, "PCA Silv" = pcaSilver))
}

.idx_hc <- function(hc, n) {
  ## list array of cluster indices at each of the n-1 merges
  idx_hc <- array(list(), c(2*n-1, 2))
  idx_hc[1:n, 1] <- as.list(n:1)
  idx_hc[(n+1):(2*n-1), ] <- hc$merge + n + (hc$merge<0)
  
  ## complete idx_hc
  for (k in 1:(n-1)) {
    idx_hc[[n+k, 1]] <- unlist(idx_hc[idx_hc[[n+k, 1]], ])
    idx_hc[[n+k, 2]] <- unlist(idx_hc[idx_hc[[n+k, 2]], ])
  }
  
  ## flip index, hclust and shc use revered ordering
  idx_hc[(2*n-1):(n+1), ]
}


#create a function to run the dip test on significant nodes returned by sigclust
.run_dip <- function(Ngenes, obj, varGenes, pcs = F){
    #sort the genes
    #var_genes_sorted <- .SortGenes(obj)
    
    #select the top N genes
    #top_N <- var_genes_sorted$gene[1:Ngenes]
    
    if(pcs == T){
      obj_sub <- obj
      obj_shc <- shc(obj_sub, matmet = .fast_dist, linkage="ward", alpha=0.05, rcpp = T, n_min = 50)
    }
    
    if(pcs == F){
      #select the top N genes
      obj_sub <- as.matrix(obj[varGenes$gene[1:Ngenes], ])
      
      #Run significant clustering
      obj_shc <- shc(t(obj_sub), matmet = .fast_dist, linkage="ward", alpha=0.05, rcpp = T, n_min = 50)
      obj_sub <- t(obj_sub)
    }

    p <- plot(obj_shc, hang=.1, main = 'obj')
    
    #compute the number of significant and non significant nodes
    sig = obj_shc$nd_type == "sig"
    sig_ind <- obj_shc$idx_hc[sig, ]
    
    #get the indices of the significant nodes
    sig_arr <- which(sig, arr.ind = TRUE, useNames = TRUE)
    not_sig <- obj_shc$nd_type == "not_sig"
    not_sig_ind <- obj_shc$idx_hc[not_sig, ]
    
    #compute the distances speed things up
    dists <- parallelDist(obj_sub)
    
    
    #perform differential expression for both significant and non significant nodes
    print("Now running dip test for significant nodes")
    #perform DE for significant nodes
    Dip_Results = NULL
    dip_stat <- NULL
    if(length(sig_ind) > 0){
      count <- 1
      x = length(sig_ind)/2
      for (i in 1:x){
        c1 <- t(obj_sub[sig_ind[[i]], ])
        c2 <- t(obj_sub[sig_ind[[i+x]], ])
        Xnew <- cbind(c1, c2)
        #print(colnames(Xnew))
        
        dist_vals <- subset(dists, colnames(Xnew))
        #print(as.matrix(dist_vals)[1:5,1:5])
        
        #run the dip test
        dip_val <- dip.test(dist_vals)$p.value
        print(dip_val)
        Dip_Results[[i]] <- dip_val
        
        #run the dip algorithm
        dip_stat[[i]] <- dip(dist_vals, full.result = T)
        
        
        
      }
    }
    return(list("Dip" = Dip_Results, "plot" = p, "Dip Stat" = dip_stat, "shc" = obj_shc))
    
}

# Create a function to compute a fast euclidean distance

.fast_dist <- function(mat){
  #compute the euclidean distance using using a parallel method
  dist_mat <- parallelDist::parallelDist(mat, method = "euclidean")
  
  #return the dist matrix
  return(dist_mat)
}


# create a function 
.perm_test <- function(df, niter= 100){
  c1 <- df[,1]
  c2 <- df[, 2]

  #for(i in 1:niter){
    c1_new <- gtools::permute(c1)
    #compute the distance between both vectors
    dist_boot <- as.numeric(parallelDist::parallelDist(rbind(c1_new, c2)))
    return(dist_boot)
 
  
}
  
.Compute_Z <- function(b, data){
  my_dist <- as.numeric(parallelDist::parallelDist(data))
  # perform a significance test on the two vectors
  # return the z score of the test
  test <- BSDA::z.test(b, mu = my_dist, sigma.x = sd(b))
  z_score <- test$statistic
  #p_val <- test$p.value
  return(z_score)
}

# Create a function to run SCRNA
.RunScrna <- function(counts, dim = 2, nGenes = 500){
  #remove all rows and columns containing all zeros
  counts<-counts[rowSums(counts>0)>0,] #remove rows with zeros all the way across
  counts<-counts[,colSums(counts>0)>0]
  
  # Get the most informative genes
  gm <- compute_gene_info(as.matrix(counts),gmeta=rownames(counts),mod="multinomial")
  
  # sort the genes
  gm_sorted <- order(gm$deviance,decreasing=TRUE,na.last=FALSE)
  
  # extract the counts containing the top n most informative genes
  counts_filtered <- counts[gm_sorted[1:nGenes], ]
  
  # run glmp-pca using the filtered dataset
  glm_res <- glmpca(as.matrix(counts_filtered), dim, fam = "poi")
  
  # compute the poisson residuals
  resids <- null_residuals(as.matrix(counts_filtered), type="deviance", mod="poisson")
  
  # run pca on the residuals
 # pca_res <- prcomp(as.matrix(t(resids)),center=F,scale.=F, rank. = dim)
  
  # return the results
  return(list("counts_filtered" = counts_filtered, "residuals"= resids, "pca_residuals" = glm_res))
}

.RunSeurat <- function(counts, nGenes = 1000, dim = 2){
  # create and preprocess the object
  obj <- CreateSeuratObject(raw.data = counts, project = "Sim", minCells = 5, minGenes = 200)
  obj <- .PreProcess(obj, low_thres = c(200, -Inf, -Inf), high_thres = c(nGenes, 0.05, Inf))
  obj <- ScaleData(obj, num.cores = 6)
  varGenes <- .SortGenes(obj)
  obj_sub <- .GetGeneData(obj, varGenes$gene)
  obj_sub <- as.data.frame(obj_sub)
  obj_sub_norm <- .GetGeneData(obj, varGenes$gene, scaled = F)
  
  # Run PCA and return the first two PCS
  obj <- RunPCA(object = obj, do.print = F, pcs.compute = dim)
  pcs <- obj@dr$pca@gene.loadings[, 1:2]
  
  # return the results
  return(list("pcs" = pcs, "obj_processed" = obj_sub, "obj_processed_norm" = obj_sub_norm))
}