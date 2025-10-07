# DE test each control gRNA against all others
de_testing_cntrl <- function(sce_cntrl, cntrl_ids) {
  
  sce_cntrl_filt <- sce_cntrl[, sce_cntrl$gRNA %in% cntrl_ids]
  
  metadata <- colData(sce_cntrl_filt) %>% 
    as.data.frame()
  
  # design
  design <- model.matrix(~ 0 + gRNA + individual + batch, data = metadata)
  is.fullrank(design)
  
  # fit LM with limma
  y <- new("EList")
  y$E <- logcounts(sce_cntrl_filt)
  fit <- limma::lmFit(object = y, design = design)
  
  # contrasts
  n <- length(cntrl_ids)
  m <- ncol(design) - n
  contrast <- sapply(0:(n - 1), function(x) {c(append(rep(-1 / (n - 1), n - 1), 1, after = x), rep(0, m))}) # append(..., 1, after = x): x chooses the gRNA to compare to the rest, this gets a contrast of 1, rep(-1 / (n - 1), n - 1): the rest of the gRNAs (we have n -1 of them) get a contrast of -1 / (n - 1), c(..., 0) : the contrast of individual is 0, sapply(0:(n-1), ...): each gRNA is chosen once to be compared to the rest
  rownames(contrast) <- colnames(design)
  colnames(contrast) <- paste0(sort(unique(metadata$gRNA)), "_VS_other")
  
  # compute coefficients for a given set of contrasts
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2, trend = T, robust = T)
  
  # get DE genes for each contrast
  de_results_cntrl <- foreach(contrast = colnames(contrast),
                              .combine = bind_rows) %do% {
                                
                                topTable(fit2, coef = contrast, number = Inf, adjust.method = "BH") %>%
                                  rownames_to_column("gene") %>% 
                                  dplyr::mutate(gRNA = gsub("_VS_other", "", contrast),
                                                n_cells = sum(metadata$gRNA == unique(gRNA)))
                                
                              }
  
  return(list(de_results = de_results_cntrl,
              fit2 = fit2))
  
}


# plot the number of DE genes per control gRNA
plot_n_de_genes <- function(de_results, p.adj_cutoff = 0.05, sigma_cutoff = 3, colors) {
  
  n_de_genes_cntrl <- de_results %>% 
    group_by(gRNA, species) %>% 
    dplyr::summarise(n_de_genes = sum(adj.P.Val < p.adj_cutoff),
                     n_cells = unique(n_cells)) %>% 
    ungroup() %>% 
    dplyr::mutate(mean_de_genes = mean(n_de_genes), sd_de_genes = sd(n_de_genes), cutoff_de_genes = mean_de_genes + sigma_cutoff*sd_de_genes,
                  n_cells2 = n_cells,
                  gRNA = gsub("hg38_macFas6_NT_", "", gRNA),
                  gRNA = factor(gRNA, paste0("n", 1:47))) %>% 
    pivot_longer(cols = c("n_cells", "n_de_genes"), names_to = "variable", values_to = "n")
  
  n_de_genes_cntrl$cutoff_de_genes[n_de_genes_cntrl$variable == "n_cells"] <- NA
  
  nudge <- max(n_de_genes_cntrl$n[n_de_genes_cntrl$variable == "n_de_genes"])*0.05 
  
  p <- ggplot(n_de_genes_cntrl, aes(x = replace_na(n, 0), y = gRNA, fill = gRNA)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, guide = "none") +
    theme_bw() +
    facet_grid(species~variable, scales = "free") +
    geom_text(aes(label = replace_na(as.character(n), "NA")), nudge_x = nudge, size = 2.5) +
    scale_y_reordered(limits = rev) +
    geom_vline(aes(xintercept = cutoff_de_genes), linetype = "dashed")+
    theme(axis.title = element_blank())
  
  return(p)
  
}


# identify control gRNAs with unexpectedly many ( > mean + 3 SD) DE genes
get_bad_cntrls <- function(de_results, p.adj_cutoff = 0.05, sigma_cutoff = 3) {
  
  de_results %>% 
    group_by(gRNA, species) %>% 
    dplyr::summarise(n_de_genes = sum(adj.P.Val < p.adj_cutoff)) %>% 
    ungroup() %>% 
    dplyr::mutate(mean_de_genes = mean(n_de_genes), sd_de_genes = sd(n_de_genes), cutoff_de_genes = mean_de_genes + sigma_cutoff*sd_de_genes) %>% 
    dplyr::filter(n_de_genes > cutoff_de_genes) %>% 
    distinct(gRNA, species)
  
}


# plot UMAP colored by various variables
plot_umap <- function(seu, color_by, colors, reduction = "umap", blend = TRUE, text_size = 30, point_size = text_size / 3000, color_filt = TRUE) {
  
  metadata <- seu@meta.data
  
  if (!"cell" %in% colnames(metadata)) {
    
    metadata <- metadata %>% 
      rownames_to_column("cell")
    
  }
  
  umap_df <- seu[[reduction]]@cell.embeddings %>% 
    as.data.frame() %>% 
    rownames_to_column("cell") %>% 
    left_join(metadata) %>% 
    dplyr::rename(umap_1 = .data[[paste0(gsub("\\.|\\_", "", reduction), "_1")]], umap_2 = .data[[paste0(gsub("\\.|\\_", "", reduction), "_2")]])
  
  if (!is.numeric(umap_df[[color_by]]) & !is.null(names(colors)) & color_filt) {
    
    colors <- colors[names(colors) %in% umap_df[[color_by]]]
    
  }
  
  if (!is.numeric(umap_df[[color_by]]) & !is.null(names(colors))) {
    
    umap_df[[color_by]] <- factor(umap_df[[color_by]], names(colors))
    
  }
  
  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = .data[[color_by]])) +
    theme_bw(base_size = text_size) +
    scale_color_manual(values = colors) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.height = unit(text_size / 25, "cm"),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  if (blend) {
    
    p <- p +
      geom_point(size = point_size, alpha = 0.3, show.legend = TRUE) %>%  partition(vars(get(color_by))) * (blend("lighten") + blend("multiply", alpha = 0.5))
    
  } else {
    
    p <- p +
      geom_point(size = point_size, alpha = 0.3, show.legend = TRUE)
    
  }
  
  if (is.numeric(umap_df[[color_by]])) {
    
    p <- p +
      scale_color_gradientn(colors = colors, name = wrap_labels(gsub("\\_", " ", color_by), 15))
    
  } else {
    
    p <- p +
      scale_color_manual(values = colors, drop = FALSE, name = wrap_labels(gsub("\\_", " ", color_by), 15)) +
      guides(color = guide_legend(override.aes = list(size = text_size / 3.75, alpha = 1)))
    
  }
  
  p
  
}

# plot UMAPs per individual/species, colored by various variables
plot_umap_split <- function(seu, color_by, split_by ="individual", colors = NULL, reduction = "umap_per_indiv", blend = TRUE, text_size = 30, point_size = text_size / 3000, legend_title = color_by) {
  
  if (color_by %in% rownames(seu)) {
    
    seu@meta.data[[color_by]] <- as(seu[["RNA"]]$data, Class = "dgCMatrix")[color_by, ]
    
    if (legend_title == color_by) legend_title <- paste0(color_by, "\nexpression\n(logcounts)")
    
  }
  
  umap_df <- seu[[reduction]]@cell.embeddings %>% 
    as.data.frame() %>% 
    rownames_to_column("cell") %>% 
    left_join(seu@meta.data) %>% 
    dplyr::rename(umap_1 = .data[[paste0(gsub("\\.|\\_", "", reduction), "_1")]], umap_2 = .data[[paste0(gsub("\\.|\\_", "", reduction), "_2")]])
  
  if (!is.numeric(umap_df[[color_by]]) & !is.null(names(colors))) {
    
    colors_filt <- colors[names(colors) %in% umap_df[[color_by]]]
    umap_df[[color_by]] <- factor(umap_df[[color_by]], names(colors_filt))
    
  }
  
  p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = .data[[color_by]])) +
    theme_bw(base_size = text_size) +
    facet_wrap(~.data[[split_by]], scales = "free") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.height = unit(text_size / 25, "cm"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing = ggplot2::unit(0.1, "lines"),
          strip.background = element_blank(),
          strip.text.x = element_blank())
  
  if (blend & !is.numeric(umap_df[[color_by]])) {
    
    p <- p +
      geom_point(size = point_size, alpha = 0.3) %>%  partition(vars(get(color_by))) * (blend("lighten") + blend("multiply", alpha = 0.5))
    
  } else {
    
    p <- p +
      geom_point(size = point_size, alpha = 0.3)
    
  }
  
  if (is.numeric(umap_df[[color_by]])) {
    
    if (is.null(colors)) colors <- c("#C1E8FF", "#003F75", "black")
    
    p <- p +
      scale_color_gradientn(colors = colors, name = legend_title)
    
  } else {
    
    if (!is.null(colors)) {
      
      p <- p +
        scale_color_manual(values = colors, name = legend_title)
      
    } else {
      
      p <- p +
        labs(color = wrap_labels(gsub("\\_", " ", color_by), 15))
      
    }
    
    p <- p +
      guides(color = guide_legend(override.aes = list(size = text_size / 3.75, alpha = 1)))
    
  }
  
  p
  
}


# plot UMAPs per individual/species, with the best POU5F1 gRNA pair highlighted
plot_umap_split2 <- function(seu, color_by, split_by ="individual", colors = NULL, reduction = "umap_per_indiv", blend = TRUE, text_size = 30, point_size = text_size / 3000, legend_title = color_by) {
  
  umap_df <- seu[[reduction]]@cell.embeddings %>% 
    as.data.frame() %>% 
    rownames_to_column("cell") %>% 
    left_join(seu@meta.data) %>% 
    dplyr::rename(umap_1 = .data[[paste0(gsub("\\.|\\_", "", reduction), "_1")]], umap_2 = .data[[paste0(gsub("\\.|\\_", "", reduction), "_2")]])
  
  gRNAs <- levels(umap_df[[color_by]])[1]
  
  ggplot(umap_df, aes(x = umap_1, y = umap_2, color = .data[[color_by]])) +
    theme_bw(base_size = text_size) +
    facet_wrap(~.data[[split_by]], scales = "free") +
    geom_point(data = . %>% dplyr::filter(.data[[color_by]] != gRNAs), size = point_size, alpha = 0.3) +
    geom_point(data = . %>% dplyr::filter(.data[[color_by]] == gRNAs), size = point_size*1.5, alpha = 0.7) +
    scale_color_manual(values = colors, breaks = c("best POU5F1\ngRNA pair", "other POU5F1\ngRNAs", "control")) +
    guides(color = guide_legend(override.aes = list(size = text_size / 3.75, alpha = 1))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.height = unit(text_size / 25, "cm"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing = ggplot2::unit(0.1, "lines"),
          legend.key.spacing.y = unit(0.5, "cm"),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank())
  
}


# downsample cells so that cell numbers are equal across species, and if possible, also across individuals and batches
downsample_cells <- function(metadata, perturbations, N = NULL) {
  
  metadata_TF <- metadata %>% 
    dplyr::filter(perturbation %in% perturbations)

  human_exp <- metadata_TF %>% 
    dplyr::filter(species == "human") %>% 
    pull(batch) %>% 
    unique()
  cynomolgus_exp <- metadata_TF %>% 
    dplyr::filter(species == "cynomolgus") %>% 
    pull(batch) %>% 
    unique()
  
  metadata_TF <- metadata_TF %>% 
    dplyr::filter(batch %in% intersect(human_exp, cynomolgus_exp)) %>% 
    dplyr::mutate(batch = factor(batch, intersect(human_exp, cynomolgus_exp)))
  
  if (is.null(N)) {
    
    N <- metadata_TF %>% 
      dplyr::count(species) %>% 
      pull(n) %>% 
      min()
    
  } 
  
  species <- setNames(levels(metadata_TF$species),
                      levels(metadata_TF$species))
  
  cell_num_per_indiv <- lapply(species, function(spec) {
    
    metadata_TF %>% 
      dplyr::filter(species == spec) %>% 
      dplyr::count(individual) %>% 
      arrange(n) %>% 
      deframe()
    
  })
  
  if (length(unlist(cell_num_per_indiv)) != 4)
    stop("One of the individuals don't have any cells.")
  
  cell_num_per_indiv_downsampl <- lapply(cell_num_per_indiv, calculate_balanced_cell_numbers, cn_total = N) %>% 
    unname() %>% 
    unlist()
  
  individuals <- setNames(unique(metadata_TF$individual),
                          unique(metadata_TF$individual))
  
  cell_num_per_gRNA <- lapply(individuals, function(indiv) {
    
    metadata_TF %>% 
      dplyr::filter(individual == indiv) %>% 
      dplyr::count(perturbation) %>% 
      arrange(n) %>% 
      deframe()
    
  })
  
  cell_num_per_gRNA_downsampl <- lapply(individuals, function(indiv) {
    
    calculate_cell_numbers_per_gRNA(cell_num_per_gRNA[[indiv]], cell_num_per_indiv_downsampl[indiv])
    
  })
  
  cell_num_per_batch <- lapply(individuals, function(indiv) {
    
    perturbations <- setNames(names(cell_num_per_gRNA_downsampl[[indiv]]),
                              names(cell_num_per_gRNA_downsampl[[indiv]]))
    
    lapply(perturbations, function(pert) {
      
      metadata_TF %>% 
        dplyr::filter(individual == indiv & perturbation == pert) %>% 
        dplyr::count(batch) %>% 
        arrange(n) %>% 
        deframe()
      
    }) 
    
  })
  
  cell_num_per_batch_downsampl <- lapply(individuals, function(indiv) {
    
    perturbations <- setNames(names(cell_num_per_gRNA_downsampl[[indiv]]),
                              names(cell_num_per_gRNA_downsampl[[indiv]]))
    
    lapply(perturbations, function(pert) {
      
      calculate_cell_numbers_per_gRNA(cell_num_per_batch[[indiv]][[pert]], cell_num_per_gRNA_downsampl[[indiv]][pert])
      
    })
    
  })
  
  cells_downsampl <- lapply(individuals, function(indiv) {
    
    cell_num_per_batch_downsampl_indiv <- cell_num_per_batch_downsampl[[indiv]]
    
    lapply(names(cell_num_per_batch_downsampl_indiv), function(pert) {
      
      cell_num_per_batch_downsampl_indiv_gRNA <- cell_num_per_batch_downsampl_indiv[[pert]]
      
      lapply(names(cell_num_per_batch_downsampl_indiv_gRNA), function(b) {
        
        cells_TF_indiv_gRNA_batch <- metadata_TF %>% 
          dplyr::filter(individual == indiv & perturbation == pert & batch == b) %>% 
          pull(cell)
        
        cn_target <- cell_num_per_batch_downsampl_indiv_gRNA[b]
        
        set.seed(0)
        
        sample(cells_TF_indiv_gRNA_batch, cn_target)
        
      })
      
    })
    
  }) %>% unlist() %>% unname()
  
  return(cells_downsampl)

}


# if possible, balance cell numbers across individuals
calculate_balanced_cell_numbers <- function (cn_indiv, cn_total) {
  
  cn_individual1 <- cn_indiv[1]
  cn_individual2 <- cn_indiv[2]
  
  if (cn_individual1 >= cn_total / 2 && cn_individual2 >= cn_total /2) {
    result = c(cn_total / 2, cn_total / 2)
  } else if (cn_individual1 + cn_individual2 > cn_total) {
    if (cn_individual1 != cn_individual2) {
      result = c(cn_individual1, cn_total -  cn_individual1) 
    }
  } else {
    result = c(cn_individual1,cn_individual2)
  }
  names(result) <- names(cn_indiv)
  result
}


# if possible, balance cell numbers across gRNAs
calculate_cell_numbers_per_gRNA <- function(cn_gRNA, cn_total) {
  
  if (sum(cn_gRNA) == cn_total) {
    
    return(cn_gRNA)
    
  } else {
    
    cn_target <- floor(cn_total / length(cn_gRNA))
    
    rare_gRNA_nums <- cn_gRNA[cn_gRNA <= cn_target]
    
    cn_remaining <- cn_total - sum(rare_gRNA_nums)
    
    abundant_gRNA_nums <- cn_gRNA[cn_gRNA > cn_target]
    
    abundant_gRNA_nums_downsampl <- split_counts_capped(cn_remaining, abundant_gRNA_nums)
    
    return(c(rare_gRNA_nums,
             abundant_gRNA_nums_downsampl))
    
  }
  
}


# make counts as even as possible, while respecting a cap (maximum available number of cells)
split_counts_capped <- function(cn_remaining, abundant_gRNA_nums) {
  
  ## Input checks
  if (!is.numeric(cn_remaining) || cn_remaining %% 1 != 0 || cn_remaining < 0)
    stop("`cn_remaining` must be a non-negative integer.")
  
  if (!is.numeric(abundant_gRNA_nums) ||
      any(abundant_gRNA_nums %% 1 != 0) ||
      any(abundant_gRNA_nums < 0))
    stop("`abundant_gRNA_nums` must be a (named) vector of non-negative integers.")
  
  n <- length(abundant_gRNA_nums)
  if (n == 0)
    stop("`abundant_gRNA_nums` must contain at least one entity.")
  
  if (sum(abundant_gRNA_nums) < cn_remaining)
    stop("Requested total exceeds the combined capacity of all entities.")
  
  ## Start with the largest even split that respects caps 
  base      <- floor(cn_remaining / n)             # ideal even share
  alloc     <- pmin(base, abundant_gRNA_nums)      # but respect each cap
  remaining <- cn_remaining - sum(alloc)           # how many cells still unplaced
  
  ## Hand out the leftover cells, 1 per round
  ## At each round we give +1 to every currently least-filled entity that
  ## still has capacity.  This guarantees the final allocations differ
  ## by at most 1 unless some caps force a wider gap.
  while (remaining > 0) {
    
    candidates <- which(alloc < abundant_gRNA_nums)  # rooms with spare capacity
    min_now    <- min(alloc[candidates])             # current lowest allocation
    least_filled <- candidates[alloc[candidates] == min_now]
    
    for (i in rev(least_filled)) {
      if (remaining == 0) break
      alloc[i]  <- alloc[i] + 1
      remaining <- remaining - 1
    }
  }
  
  names(alloc) <- names(abundant_gRNA_nums)
  alloc
}


# filter and scran-normalize SCE object
filter_and_normalize <- function(sce) {
  
  # add perturbation_species concatenated variable
  sce$perturbation_species <- paste0(sce$perturbed_TF, "_", sce$species)
  
  # get (subsetted) metadata and raw counts
  metadata <- as.data.frame(colData(sce))
  cnts <- counts(sce)
  
  # get expressed genes
  gene_list <- foreach(id = unique(sce$perturbation_species)) %do% {
    
    indivs <- metadata %>%
      dplyr::filter(perturbation_species == id) %>%
      pull(individual) %>%
      unique()
    
    cnt_mat_list_id <- lapply(indivs, function(indiv) {
      
      cnts[,sce$perturbation_species == id & sce$individual == indiv]
      
    })
    
    gene_list_id <- lapply(cnt_mat_list_id, function(cnt_mat) {
      
      gene_logical <- (rowSums(cnt_mat > 0) / ncol(cnt_mat)) >= 0.1 & rowSums(cnt_mat > 0) >= 10
      names(gene_logical)[gene_logical]
      
    })
    
    Reduce(intersect, gene_list_id)
    
  }
  expr_genes <- Reduce(union, gene_list)
  sce <- sce[expr_genes, ]
  
  # normalize
  sce <- computeSumFactors(sce)
  preclusters <- quickCluster(sce, min.size = 50)
  sce <- computeSumFactors(sce, clusters = preclusters)
  sce <- scater::logNormCounts(sce)
  
  return(sce)
  
}


# perform DE testing between conditions and species using dream (mixed effects modeling)
DE_testing_LMM <- function(sce, n_cores = 20) {
  
  # design
  if (length(unique(sce$batch)) > 1) {
    
    form <- ~ species + condition + species:condition + (1|individual) + (1|batch)
    
  } else {
    
    form <- ~ species + condition + species:condition + (1|individual)
    
  }
  
  # metadata
  metadata <- as.data.frame(colData(sce))
  TF <- setdiff(unique(metadata$perturbed_TF), "NT_control")
  
  # contrast
  L <- makeContrastsDream(form, 
                          metadata,
                          contrasts = c(DE_cyno = paste0("condition", TF, " + speciescynomolgus:condition", TF))) 
  
  # fit the dream model on each gene
  fit <- dream(as.matrix(logcounts(sce)), form, metadata, L, useWeights = FALSE, BPPARAM = MulticoreParam(n_cores), computeResiduals = FALSE)
  fit <- eBayes(fit)
  
  # get DE genes for each contrast
  de_results <- bind_rows(human = limma::topTable(fit, coef = paste0("condition", TF), number = Inf, adjust.method = "BH", confint = T) %>% 
                            rownames_to_column("gene"),
                          cynomolgus = limma::topTable(fit, coef = "DE_cyno", number = Inf, adjust.method = "BH", confint = T) %>%
                            rownames_to_column("gene"),
                          interaction = limma::topTable(fit, coef = paste0("speciescynomolgus:condition", TF),  number = Inf, adjust.method = "BH", confint = T) %>% 
                            rownames_to_column("gene"),
                          .id = "contrast")
  
  return(de_results)
  
}


# plot DE results across positively or negatively correlated module memebr genes and non-module genes
plot_DE_results <- function(de_results, pruned_modules, all_genes) {
  
  # all genes tested in Perturb-seq DE analysis
  all_genes_perturbseq <- unique(de_results$gene)
  
  # intersection of genes in the network and genes in Perturb-seq data
  shared_genes <- intersect(all_genes_perturbseq, all_genes)
  
  # filter objects
  pruned_modules_filt <- pruned_modules %>% 
    dplyr::filter(target %in% shared_genes & regulator == "POU5F1") %>% 
    dplyr::select(regulator, target, category = direction)
  
  de_results_filt <- de_results %>% 
    dplyr::filter(gene %in% shared_genes & gene != "POU5F1") %>% 
    left_join(pruned_modules_filt,
              by = c("gene" = "target")) %>% 
    dplyr::mutate(category = case_when(category == "+" ~ "activated\ntargets",
                                       category == "-" ~ "repressed\ntargets",
                                       is.na(category) ~ "not in\nmodule")) %>% 
    dplyr::mutate(category = factor(category, levels = c("not in\nmodule", "activated\ntargets", "repressed\ntargets"))) %>% 
    dplyr::select(category, contrast, gene, logFC, p_adj = adj.P.Val, logFC_lwr = CI.L, logFC_upr = CI.R) %>%
    pivot_wider(names_from = "contrast", values_from = c("logFC", "logFC_lwr", "logFC_upr", "p_adj"), names_sep = ".")
  
  # plot
  max_lfc <- max(c(de_results_filt$logFC_upr.human, de_results_filt$logFC_upr.cynomolgus))
  min_lfc <- min(c(de_results_filt$logFC_lwr.human, de_results_filt$logFC_lwr.cynomolgus))
  
  de_results_filt %>%
    ggplot(aes(x = logFC.human, y = logFC.cynomolgus, color = category, size = category)) +
    geom_errorbar(aes(xmin = logFC_lwr.human, xmax = logFC_upr.human, color = category), size = 0.07, alpha = 0.5, show.legend = FALSE) +
    geom_errorbar(aes(ymin = logFC_lwr.cynomolgus, ymax = logFC_upr.cynomolgus, color = category), size = 0.07, alpha = 0.5, show.legend = FALSE) +
    geom_point(data = . %>% dplyr::filter(category == "not in\nmodule"), aes(color = category), size = 0.005, alpha = 0.7) +
    geom_point(data = . %>% dplyr::filter(category != "not in\nmodule"), aes(color = category), size = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    theme_bw(base_size = 14) +
    scale_color_manual(values = c("activated\ntargets" = "darkseagreen4", "repressed\ntargets" = "indianred3", "not in\nmodule" = "grey70")) +
    scale_size_manual(values = c("TRUE" = 0.2, `no significant difference` = 0.05), guide = "none") +
    geom_label_repel(data = . %>% filter(gene == "SPP1"),
                     aes(label = gene), fill = "white", size = 3, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, max.overlaps = 10, show.legend = FALSE) +
    xlab(expression(paste(log[2], FC, " (human)"))) +
    ylab(expression(paste(log[2], FC, " (cynomolgus)"))) +
    theme(legend.title = element_blank(),
          legend.key.spacing.y = unit(0.4, "cm"),
          legend.margin = margin(l = -10)) +
    xlim(min_lfc, max_lfc) +
    ylim(min_lfc, max_lfc)
  
}


# plot DE and DR results and highlight module member genes
plot_module_members <- function(de_results, pruned_modules, all_genes) {
  
  # all genes tested in Perturb-seq DE analysis
  all_genes_perturbseq <- unique(de_results$gene)
  
  # intersection of genes in the network and genes in Perturb-seq data
  shared_genes <- intersect(all_genes_perturbseq, all_genes)
  
  # filter objects
  pruned_modules_filt <- pruned_modules %>% 
    dplyr::filter(target %in% shared_genes & regulator == "POU5F1") %>% 
    dplyr::select(regulator, target, category = direction)
  
  de_results_filt <- de_results %>% 
    dplyr::filter(gene %in% shared_genes & contrast != "interaction" & gene != "POU5F1") %>% 
    left_join(pruned_modules_filt,
              by = c("gene" = "target")) %>% 
    dplyr::mutate(category = case_when(category == "+" ~ "activated\ntargets",
                                       category == "-" ~ "repressed\ntargets",
                                       is.na(category) ~ "not in\nmodule")) %>% 
    dplyr::mutate(contrast = factor(contrast, levels = c("human", "cynomolgus")),
                  category = factor(category, levels = c("not in\nmodule", "activated\ntargets", "repressed\ntargets")))
  
  # POU5F1
  de_results_filt %>% 
    dplyr::mutate(x = paste0(contrast, "_", category),
                  x = factor(x, levels = c("human_not in\nmodule", "human_activated\ntargets", "human_repressed\ntargets", "spacer", "cynomolgus_not in\nmodule", "cynomolgus_activated\ntargets", "cynomolgus_repressed\ntargets"))) %>% 
    ggplot(aes(x = x, y = logFC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_boxplot(aes(fill = category), linewidth = 0.2, outlier.size = 0.1, outlier.alpha = 0.7) +
    theme_bw(base_size = 14) +
    scale_fill_manual(values = c("activated\ntargets" = "darkseagreen4", "repressed\ntargets" = "indianred3", "not in\nmodule" = "grey70")) +
    ylab(expression(paste(log[2], "FC")))  +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 14, color = "black"),
          legend.title = element_blank()) +
    geom_signif(comparisons = list(c("human_not in\nmodule", "human_activated\ntargets"), 
                                   c("human_not in\nmodule", "human_repressed\ntargets"),
                                   c("cynomolgus_not in\nmodule", "cynomolgus_activated\ntargets"), 
                                   c("cynomolgus_not in\nmodule", "cynomolgus_repressed\ntargets")),
                map_signif_level = c("****"=0.0002, "***"=0.002, "**"=0.02, "*"=0.1, "." = 0.2), # ugly implementation of one-tailed tests
                y_position = c(3.2, 3.45, 3.2, 3.45),
                textsize = 2.5, size = 0.2, vjust = c(rep(0.4, 3), rep(0.4, 3), rep(0.4, 6)), tip_length = 0.02) +
    scale_x_discrete(limits =  c("human_not in\nmodule", "human_activated\ntargets", "human_repressed\ntargets", "spacer", "cynomolgus_not in\nmodule", "cynomolgus_activated\ntargets", "cynomolgus_repressed\ntargets"), breaks = c("human_activated\ntargets", "cynomolgus_activated\ntargets"), labels = c("human", "cynomolgus"))
  
}


# plot POU5F1 and SPP1 expression levels
plot_POU5F1_SPP1_expr <- function(sce) {
  
  logcnts <- logcounts(sce)
  
  sce %>% 
    colData() %>% 
    as.data.frame() %>% 
    dplyr::mutate(condition = factor(ifelse(perturbed_TF == "NT_control", "control", "perturbed"), c("perturbed", "control")),
                  POU5F1 = logcnts["POU5F1",],
                  SPP1 = logcnts["SPP1",]) %>% 
    pivot_longer(cols = c("POU5F1", "SPP1"), names_to = "gene", values_to = "expr") %>% 
    ggplot(aes(y = condition, x = expr, fill = condition)) +
    geom_boxplot(linewidth = 0.2, outlier.size = 0.1, outlier.alpha = 0.7) +
    theme_bw(base_size = 14) +
    facet_nested(gene + species ~ ., scales = "free_y") +
    xlab("expression [logcounts]") +
    scale_fill_manual(values = c("orange2", "peachpuff1"), guide = "none") +
    theme(axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_blank(),
          panel.spacing.y = unit(0.01, "npc"),
          strip.text = element_text(size = 9))
  
  
}
