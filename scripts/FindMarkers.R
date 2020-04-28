.findMarkers <- function(x, clusters, design=NULL, pval.type=c("any", "all"), direction=c("any", "up", "down"), subset.row=NULL)
  # Uses limma to find the markers that are differentially expressed between clusters,
  # given a log-expression matrix and some blocking factors in 'design'.
  #
  # written by Aaron Lun
  # created 22 March 2017
  # last modified 4 May 2017    
{
  # Creating a design matrix.
  clusters <- as.factor(clusters)
  full.design <- model.matrix(~0 + clusters)
  colnames(full.design) <- clust.vals <- levels(clusters)
  
  if (!is.null(design)) {
    # Removing terms to avoid linearly dependencies on the intercept.
    out <- qr.solve(design, cbind(rep(1, nrow(design))))
    to.drop <- abs(out) > 1e-8
    if (any(to.drop)) {
      design <- design[,-which(to.drop)[1],drop=FALSE]
    }
    full.design <- cbind(full.design, design) # Other linear dependencies will trigger warnings.
  }
  
  pval.type <- match.arg(pval.type) 
  direction <- match.arg(direction)
  subset.row <- .subset_to_index(subset.row, x, byrow=TRUE)
  lfit <- limma::lmFit(x[subset.row,,drop=FALSE], full.design)
  output <- vector("list", length(clust.vals))
  names(output) <- clust.vals  
  
  for (host in clust.vals) { 
    not.host <- clust.vals!=host
    targets <- clust.vals[not.host]
    all.p <- all.lfc <- vector("list", length(targets))
    names(all.p) <- names(all.lfc) <- targets
    
    con <- matrix(0, ncol(full.design), length(clust.vals))
    diag(con) <- -1
    con[which(!not.host),] <- 1
    con <- con[,not.host,drop=FALSE]
    colnames(con) <- targets
    
    fit2 <- contrasts.fit(lfit, con)
    fit2 <- eBayes(fit2, trend=TRUE, robust=TRUE)
    
    for (target in targets) { 
      res <- topTable(fit2, number=Inf, sort.by="none", coef=target)
      pvals <- res$P.Value
      
      if (direction=="up") {
        pvals <- ifelse(res$logFC > 0, pvals/2, 1-pvals/2)                
      } else if (direction=="down") {
        pvals <- ifelse(res$logFC < 0, pvals/2, 1-pvals/2)                
      }
      
      all.p[[target]] <- pvals
      all.lfc[[target]] <- res$logFC
    }
    
    com.p <- do.call(rbind, all.p)
    ngenes <- ncol(com.p)
    if (pval.type=="any") { 
      # Computing Simes' p-value in a fully vectorised manner.
      ncon <- nrow(com.p)
      gene.id <- rep(seq_len(ngenes), each=ncon)
      penalty <- rep(ncon/seq_len(ncon), ngenes) 
      o <- order(gene.id, com.p)
      com.p[] <- com.p[o]*penalty
      com.p <- t(com.p)
      smallest <- (max.col(-com.p) - 1) * ngenes + seq_len(ngenes)
      pval <- com.p[smallest]
    } else {
      # Computing the IUT p-value.
      com.p <- t(com.p)
      largest <- (max.col(com.p) - 1) * ngenes + seq_len(ngenes)
      pval <- com.p[largest]
    }
    
    collected.ranks <- lapply(all.p, rank, ties="first")
    min.rank <- do.call(pmin, collected.ranks)
    names(all.lfc) <- paste0("logFC.", names(all.lfc))
    marker.set <- data.frame(Top=min.rank, Gene=rownames(x)[subset.row], 
                             FDR=p.adjust(pval, method="BH"), do.call(cbind, all.lfc), 
                             stringsAsFactors=FALSE, check.names=FALSE)
    marker.set <- marker.set[order(marker.set$Top),]
    rownames(marker.set) <- NULL
    output[[host]] <- marker.set
  }
  
  return(output)
}

setGeneric("findMarkers", function(x, ...) standardGeneric("findMarkers"))

setMethod("findMarkers", "matrix", .findMarkers)

setMethod("findMarkers", "SCESet", function(x, ..., subset.row=NULL, assay="exprs", get.spikes=FALSE) {
  if (is.null(subset.row)) { subset.row <- .spike_subset(x, get.spikes) }
  .findMarkers(assayDataElement(x, assay), ..., subset.row=subset.row)
})  


.subset_to_index <- function(subset, x, byrow=TRUE) {
  if (byrow) {
    dummy <- seq_len(nrow(x))
    names(dummy) <- rownames(x)
  } else {
    dummy <- seq_len(ncol(x))
    names(dummy) <- colnames(x) 
  }
  
  if (!is.null(subset)) { 
    dummy <- dummy[subset]
  }
  out <- unname(dummy)
  if (any(is.na(out))) {
    stop("'subset' indices out of range of 'x'")
  }
  return(out)
}



get_auroc <- function(gene, labels) {
  score <- rank(gene)
  # Get average score for each cluster
  ms <- aggregate(score ~ labels, FUN = mean)
  # Get cluster with highest average score
  posgroup <- ms[ms$score == max(ms$score), ]$labels
  # Return NAs if there is a tie for cluster with highest average score (by definition this is
  # not cluster specific)
  if (length(posgroup) > 1) {
    return(c(NA, NA, NA))
  }
  # Create 1/0 vector of truths for predictions, cluster with highest average score vs
  # everything else
  truth <- as.numeric(labels == posgroup)
  # Make predictions & get auc using RCOR package.
  pred <- ROCR::prediction(score, truth)
  val <- unlist(ROCR::performance(pred, "auc")@y.values)
  pval <- suppressWarnings(wilcox.test(score[truth == 1], score[truth == 0])$p.value)
  return(c(val, posgroup, pval))
}



#' Calculate marker genes
#'
#' Find marker genes in the dataset. The \code{\link{get_auroc}} is used to calculate
#' marker values for each gene.
#'
#' @param dataset expression matrix
#' @param labels cell labels corresponding clusters
#' @return data.frame containing the marker genes, corresponding cluster indexes
#' and adjusted \code{p-value}s
#' @importFrom stats p.adjust
#' @examples
#' d <- get_marker_genes(yan[1:10,], as.numeric(ann[,1]))
#' d
#' 
#' @export
get_marker_genes <- function(dataset, labels) {
  res <- apply(dataset, 1, get_auroc, labels = labels)
  res <- data.frame(matrix(unlist(res), ncol = 3, byrow = T))
  colnames(res) <- c("auroc", "clusts", "pvalue")
  rownames(res) <- rownames(dataset)
  res$pvalue <- p.adjust(res$pvalue)
  return(res)
}


#' Reorder and subset gene markers for plotting on a heatmap
#' 
#' Reorders the rows of the input data.frame based on the \code{sc3_k_markers_clusts}
#' column and also keeps only the top 10 genes for each value of \code{sc3_k_markers_clusts}.
#'
#' @param markers a \code{data.frame} object with the following colnames:
#' \code{sc3_k_markers_clusts}, \code{sc3_k_markers_auroc}, \code{sc3_k_markers_padj}.
#' 
markers_for_heatmap <- function(markers) {
  res <- NULL
  for (i in unique(markers[, 2])) {
    tmp <- markers[markers[, 2] == i, ]
    if (nrow(tmp) > 10) {
      res <- rbind(res, tmp[1:10, ])
    } else {
      res <- rbind(res, tmp)
    }
  }
  
  return(res)
}


#' Plot expression of marker genes identified by \code{SC3} as a heatmap.
#' 
#' By default the genes with the area under the ROC curve (AUROC) > 0.85 
#' and with the p-value < 0.01 are selected and the top 10 marker 
#' genes of each cluster are visualized in this heatmap.
#' 
#' @name sc3_plot_markers
#' @aliases sc3_plot_markers, sc3_plot_markers,SingleCellExperiment-method
#' 
#' @param object an object of 'SingleCellExperiment' class
#' @param k number of clusters
#' @param auroc area under the ROC curve
#' @param p.val significance threshold used for the DE genes
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#' 
#' @importFrom pheatmap pheatmap
Plot_Markers <- function(dat, auroc = 0.85, p_val = 0.05, labels ) {
  
  add_ann_col <- FALSE
  ann <- NULL
  # if (!is.null(show_pdata)) {
  #   ann <- make_col_ann_for_heatmaps(object, show_pdata)
  #   if (!is.null(ann)) {
  #     add_ann_col <- TRUE
  #     # make same names for the annotation table
  #     rownames(ann) <- colnames(dataset)
  #   }
  # }
  
  # get all marker genes
  k <- length(unique(labels))
  markers <- get_marker_genes(dat, ensemble.labs)
  markers <- organise_marker_genes(markers, k, p_val, auroc)
  
  if(!is.null(markers)) {
    # get top 10 marker genes of each cluster
    markers <- markers_for_heatmap(markers)
    
    row.ann <- data.frame(Cluster = factor(markers[, 2], levels = unique(markers[, 2])))
    rownames(row.ann) <- rownames(markers)
    
    do.call(pheatmap::pheatmap, c(list(dat[rownames(markers), , drop = FALSE], show_colnames = FALSE, 
                                       cluster_rows = FALSE, cluster_cols = labels, cutree_cols = k, annotation_row = row.ann, annotation_names_row = FALSE, 
                                       gaps_row = which(diff(markers[, 1]) != 0), cellheight = 10), list(annotation_col = ann)[add_ann_col]))
  } else {
    message("No markers have been found, try to lower significance thresholds!")
  }
}


organise_marker_genes <- function(dat, k, p_val = 0.05, auroc = 0.85) {
  
  dat <- dat[dat[, "pvalue"] < p_val, ]
  dat <- dat[dat[, "auroc"] > auroc, ]
  
  d <- NULL
  
  for (i in sort(unique(dat[, "clusts"]))) {
    tmp <- dat[dat[, "clusts"] == i, ]
    tmp <- tmp[order(tmp[, "auroc"], decreasing = TRUE), ]
    d <- rbind(d, tmp)
  }
  
  if(nrow(dat) > 0) {
    return(d)
  } else {
    return(NULL)
  }
}


Plot.DE.Genes <- function(dat, genes, p.val = 0.05, labels) {
  
  
  add_ann_col <- FALSE
  ann <- NULL
  # if (!is.null(show_pdata)) {
  #   ann <- make_col_ann_for_heatmaps(object, show_pdata)
  #   if (!is.null(ann)) {
  #     add_ann_col <- TRUE
  #     # make same names for the annotation table
  #     rownames(ann) <- colnames(dataset)
  #   }
  # }
  
  # de_genes <- organise_de_genes(object, k, p.val)
  de_genes <- head(genes, 50)
  # remove Inf when the p-value is actually 0 (less than the accuracy limit)
  de_genes[de_genes < 1e-17] <- 1e-17
  row_ann <- data.frame(log10_padj = -log10(de_genes))
  rownames(row_ann) <- names(de_genes)
  k <- length(unique(labels))
  do.call(pheatmap::pheatmap, c(list(dat[names(de_genes), , drop = FALSE], show_colnames = FALSE, 
                                     cluster_rows = FALSE, cluster_cols = labels, cutree_cols = k, annotation_row = row_ann, cellheight = 10), 
                                list(annotation_col = ann)[add_ann_col]))
}


get_de_genes <- function(dataset, labels) {
  tmp <- apply(dataset, 1, kruskal.test, g = factor(labels))
  ps <- unlist(lapply(tmp, "[[", "p.value"))
  ps <- p.adjust(ps)
  return(ps)
}

