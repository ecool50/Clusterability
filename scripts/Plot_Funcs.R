##pull label information for when hang > 0
## necessary since ggdendro package won't provide this output
.lab_height <- function(tree, heights = c()) {
  for (k in seq(length(tree))) {
    if (is.leaf(tree[[k]])) {
      heights <- c(heights, attr(tree[[k]], "height"))
    } else {
      heights <- .lab_height(tree[[k]], heights)
    }
  }
  heights
}


.PlotRes <- function(hc, type = "All", alpha = 0.001, use_labs = F, datName){
  
  
  if(type == "All"){
    ##using ggdendro package
    shc_dend <- as.dendrogram(hc$hc_dat, hang=0.5)
    shc_dendat <- ggdendro::dendro_data(shc_dend)
    shc_segs <- ggdendro::segment(shc_dendat)
    shc_labs <- ggdendro::label(shc_dendat)
    
    ##change label hight to be correct
    shc_labs$y <- .lab_height(shc_dend)
    
    fwer <- T
    # if (!is.null(groups) & length(groups) == n) {
    #   shc_labs$clusters <- groups[shc$hc_dat$order]
    # }
    
    ## flip index, hclust/dend and shc use revered ordering
    rev_nd_type <- rev(hc$nd_type)
    rev_p_use <- rev(hc$p_emp)
    
    ## significant nodes
    sig_linkvals <- as.factor(hc$hc_dat$height[rev_nd_type == "sig"])
    sig_segs <- filter(shc_segs, as.factor(y) %in% sig_linkvals)
    
    ## non-signficant nodes
    notsig_linkvals <- as.factor(hc$hc_dat$height[rev_nd_type == "not_sig"])
    notsig_segs <- filter(shc_segs, as.factor(y) %in% notsig_linkvals)
    
    ## tested nodes
    test_linkvals <- as.factor(hc$hc_dat$height[which(rev_nd_type == "sig")])
    test_segs <- filter(shc_segs, as.factor(y) %in% test_linkvals)
    test_segtops <- filter(test_segs, (y == yend) & (x < xend))
    test_segtops <- test_segtops[order(test_segtops[, 2]), ]
    test_segtops <- cbind(test_segtops, 
                          "pval"=format(rev_p_use[which(rev_nd_type == "sig")], 
                                        digits=3, scientific=TRUE))
    
    ## small sample nodes
    skip_linkvals <- as.factor(hc$hc_dat$height[rev_nd_type == "n_small"])
    skip_segs <- filter(shc_segs, as.factor(y) %in% skip_linkvals)
    
    ## un-tested nodes (non-significant parent node)
    if (fwer) {
      fwer_linkvals <- as.factor(hc$hc_dat$height[rev_nd_type == "no_test"])
      fwer_segs <- filter(shc_segs, as.factor(y) %in% fwer_linkvals)
    }
    
    ## calculate various plotting dimensions prior to actual ggplot call
    ax_x_ref <- max(hc$hc_dat$height)
    
  } 
  
  ax_x_top <- max(ax_x_ref)*1.25
  ax_x_bot <- -max(ax_x_ref)/4
  ax_x_scale <- floor(log10(ax_x_top))
  ax_y_range <- max(shc_segs$y) - min(shc_segs$y)  
  
  
  ## make initial ggdendro with null colo
  ## colors to be used in figure
  col_tx_sig <- "#A60A3C"
  col_nd_null <- "gray"
  plot_dend <- ggplot() + 
    geom_segment(data=shc_segs, 
                 aes(x=x, y=y, xend=xend, yend=yend), 
                 color=col_nd_null) +
    theme_bw()
  
  ## add appropriate title
  plot_dend <- plot_dend + 
    ggtitle(paste(datName,": Showing all p-values below", 
                  alpha, "cutoff")) + theme(plot.title = element_text(hjust = 0.5, size=15), axis.title=element_text(size=14))
  
  
  ##add text labels for each object if desired
  if (use_labs) {
    shc_labs$txt_y <- shc_labs$y - (2-is.null(groups))*ax_y_range/30
    plot_dend <- plot_dend +
      geom_text(data=shc_labs, 
                aes(x=x, y=txt_y, label=label, 
                    hjust=1, vjust=.5, angle=90), 
                size=3)
  }
  
  ##add colored segments to plot if any branches were significant
  if (sum(rev_nd_type == "sig") > 0) {
    plot_dend <- plot_dend +
      geom_segment(data=sig_segs,
                   aes(x=x, y=y, xend=xend, yend=yend, color="sig"), 
                   size=1)
  }
  
  
  ##add p-values for segments if they were tested
  if (sum(rev_nd_type == "sig") > 0) {
    plot_dend <- plot_dend +
      geom_text(data=test_segtops, 
                aes(x=x, y=y, label=pval, hjust=-0.2, vjust=-0.5),
                col=col_tx_sig, size=4)
  }
  
  
  ##add FWER controlled segments if any branches were controlled
  
  if (fwer & sum(rev_nd_type == "no_test") > 0) {
    plot_dend <- plot_dend +
      geom_segment(data=fwer_segs,
                   aes(x=x, y=y, xend=xend, yend=yend, color="no_test"), 
                   size=1)
  }
  
  ##add colored segments if any branches had size < n_min
  if (sum(rev_nd_type == "n_small") > 0) {
    plot_dend <- plot_dend +
      geom_segment(data=skip_segs,
                   aes(x=x, y=y, xend=xend, yend=yend, color="n_small"), 
                   size=1)
  }
  
  ## add colored segments if any branches had non-significant calls
  if (sum(rev_nd_type == "not_sig") > 0) {
    plot_dend <- plot_dend +
      geom_segment(data=notsig_segs,
                   aes(x=x, y=y, xend=xend, yend=yend, color="not_sig"), 
                   size=1)
  }
  
  plot_dend <- plot_dend +
    scale_y_continuous(name="linkage", expand=c(.25, 0),
                       breaks=seq(0, ax_x_top, by=10^ax_x_scale)) +
    scale_x_continuous(name="", expand=c(.1, 0),
                       breaks=c(),
                       labels=c())
  
  return(plot_dend)
}


.MultiCutTRee <- function(obj, alpha = 0.05, ci_idx = 1, ci_emp = T, nmin = 10) {
  
  
  ## for easier calling, subset and reverse order
  if (ci_emp) {
    p_use <- obj$p_emp
  } else {
    p_use <- obj$p_norm[, ci_idx]
  }
  
  n <- length(p_use)  + 1
  
  ## determine significant clusters
  # cutoff <- .fwer_cutoff.matrix(obj$idx_hc, alpha)
  # pd_map <- .pd_map(obj$hc_dat, n)
  
  # nd_type <- rep("", n-1)
  # for (k in 1:(n-1)) {
  #   ## check if subtree is large enough
  #   if (length(unlist(obj$idx_hc[k, ])) < nmin) {
  #     nd_type[k] <- "n_small"
  #     next
  #   }
  #   ## check if parent was significant
  #   if ((k > 1) && (nd_type[pd_map[k]] != "sig")) {
  #     nd_type[k] <- "no_test"
  #     next
  #   }
  #   ## compare against p-value cutoff
  #   if (alpha < 1) {
  #     nd_type[k] <- ifelse(p_use[k] < 0.05, "sig", "not_sig")
  #   }
  # }
  
  ## determine non-sig nodes with sig parent nodes
  cl_idx <- which((obj$pd_map %in% which(obj$nd_type == "sig")) & (obj$nd_type != "sig"))
  if (length(cl_idx) == 0){
    cl_idx <- 1
    }
  
  print(cl_idx)
  clusters <- lapply(apply(obj$idx_hc[cl_idx, , drop=FALSE], 1, c), unlist)
  
  ## verify that the clusters contain all observations
  if (length(unique(unlist(clusters))) != n) {
    stop("length of cluster labels != n")
  }

  ## construct cluster label vector
  labs <- rep(-1, n)
  for (i in 1:length(clusters)) {
    labs[clusters[[i]]] <- i
  }
  
  return(labs)
}


