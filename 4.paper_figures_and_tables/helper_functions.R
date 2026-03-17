plotPruningExample <- function(moduleOfInterest, initial_modules, consensus_network, n_iterations) {
  
  # adjacency matrix
  adjMat <- as_adjacency_matrix(consensus_network, attr = "weight")
  
  # filter module data
  initial <- initial_modules %>%
    dplyr::filter(regulator == moduleOfInterest)
  
  #pruning nr. 1
  pruned1 <- initial %>%
    arrange(desc(weight)) %>%
    dplyr::mutate(rank_IS = 1:length(target),
                  cumsum_IS = cumsum(weight),
                  rank_cut_IS = uik(rank_IS, cumsum_IS),
                  IS_cut = weight[rank_IS == rank_cut_IS]) %>%
    dplyr::mutate(is_in_final_module = rank_IS <= rank_cut_IS)
  
  pruning1<- ggplot(pruned1, aes(x = rank_IS, y = cumsum_IS)) +
    geom_line() +
    geom_vline(aes(xintercept = rank_cut_IS), linetype = "dashed", linewidth = 0.5) +
    geom_area(aes(fill = is_in_final_module), alpha = 0.3) +
    scale_fill_manual(values = c("#BDEDDB", "#E6FFF6"), limits = c(T, F), guide = "none") +
    xlab(expression(paste("target rank based on ", italic(adj)["regulator"]))) +
    ylab(expression(paste("cumsum(", italic(adj)["regulator"], ")"))) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(size = 13, margin = margin(b = 1, l = -10)),
          plot.subtitle = element_text(margin = margin(0, 5.5, 5.5, -10)),
          plot.margin = margin(5.5, 46, 5.5, 5.5)) +
    labs(title = expression(paste(bold("Initial POU5F1 module "), "(", italic(M)["initial"], ")")),
         subtitle = paste0("number of genes: ", nrow(pruned1)))
  
  # pruning nr.2
  pruned2 <- pruned1 %>%
    dplyr::filter(is_in_final_module) %>%
    dplyr::mutate(kIM = Matrix::rowSums(adjMat[target, c(target, unique(as.character(regulator)))])) %>%
    arrange(desc(kIM)) %>%
    dplyr::mutate(rank_kIM = 1:length(target),
                  cumsum_kIM = cumsum(kIM),
                  rank_cut_kIM = uik(rank_kIM, cumsum_kIM),
                  kIM_cut = kIM[rank_kIM == rank_cut_kIM]) %>%
    dplyr::mutate(is_in_final_module = rank_kIM <= rank_cut_kIM)
  
  pruning2 <- ggplot(pruned2, aes(x = rank_kIM, y = cumsum_kIM)) +
    geom_line() +
    geom_vline(aes(xintercept = rank_cut_kIM), linetype = "dashed", linewidth = 0.5) +
    geom_area(aes(fill = is_in_final_module), alpha = 0.3) +
    scale_fill_manual(values = c("#91D3BC", "#BDEDDB"), limits = c(T, F), guide = "none") +
    xlab(expression(paste("target rank based on ", italic(kIM)["pruned1"]))) +
    ylab(expression(paste("cumsum(", italic(kIM)["pruned1"], ")"))) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(size = 13, margin = margin(b = 0, l = -10)),
          plot.subtitle = element_text(margin = margin(0, 5.5, 5.5, -10)),
          plot.margin = margin(5.5, 46, 5.5, 5.5)) +
    labs(title = expression(paste(bold("Intermediate module after 1st pruning step "), "(", italic(M)["pruned1"], ")")),
         subtitle = paste0("number of genes: ", nrow(pruned2)))
  
  # pruning nr. 3
  pruned3 <- pruned2 %>%
    dplyr::filter(is_in_final_module) %>%
    arrange(desc(weight)) %>%
    dplyr::mutate(rank_IS = 1:length(target),
                  cumsum_IS = cumsum(weight),
                  rank_cut_IS = uik(rank_IS, cumsum_IS),
                  IS_cut = weight[rank_IS == rank_cut_IS]) %>%
    dplyr::mutate(is_in_final_module = rank_IS <= rank_cut_IS)
  
  pruning3 <- ggplot(pruned3, aes(x = rank_IS, y = cumsum_IS)) +
    geom_line() +
    geom_vline(aes(xintercept = rank_cut_IS), linetype = "dashed", linewidth = 0.5) +
    geom_area(aes(fill = is_in_final_module), alpha = 0.4) +
    scale_fill_manual(values = c("#60AA98", "transparent"), limits = c(T, F), guide = "none") +
    new_scale_fill() +
    geom_area(aes(fill = weight > IS_cut, alpha = weight > IS_cut)) +
    scale_fill_manual(values = c("transparent", "#91D3BC"), limits = c(T, F), guide = "none") +
    scale_alpha_manual(values = c(0, 0.3), limits = c(T, F), guide = "none") +
    xlab(expression(paste("target rank based on ", italic(adj)["regulator"]))) +
    ylab(expression(paste("cumsum(", italic(adj)["regulator"], ")"))) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(size = 13, margin = margin(b = 0, l = -10)),
          plot.subtitle = element_text(margin = margin(0, 5.5, 5.5, -10)),
          plot.margin = margin(5.5, 46, 5.5, 5.5)) +
    labs(title = expression(paste(bold("Intermediate module after 2nd pruning step "), "(", italic(M)["pruned2"], ")")),
         subtitle = paste0("number of genes: ", nrow(pruned3)))
  
  if (n_iterations == 3) {
    
    final <- pruned3 %>%
      dplyr::filter(is_in_final_module)
    
    pruning4   <- ggplot() +
      labs(title = expression(paste(bold("Final pruned module "), "(", italic(M)["final"], ")")),
           subtitle = paste0("number of genes: ", nrow(final))) +
      theme_bw(base_size = 13) +
      theme(panel.background = element_blank(),
            plot.margin = margin(5.5, 46, -500, 5.5),
            plot.title = element_text(size = 13, margin = margin(b = 1, l = -10)),
            plot.subtitle = element_text(margin = margin(0, 5.5, 5.5, -10)))
    
    return(pruning1  / plot_spacer() / pruning2  / plot_spacer() / pruning3  / plot_spacer() / pruning4 + plot_layout(heights = c(1, .08, 1, .08, 1, .08, 0.3)))
    
  } else {
    
    # pruning nr.4
    pruned4 <- pruned3 %>%
      dplyr::filter(is_in_final_module) %>%
      dplyr::mutate(kIM = Matrix::rowSums(adjMat[target, c(target, unique(as.character(regulator)))])) %>%
      arrange(desc(kIM)) %>%
      dplyr::mutate(rank_kIM = 1:length(target),
                    cumsum_kIM = cumsum(kIM),
                    rank_cut_kIM = uik(rank_kIM, cumsum_kIM),
                    kIM_cut = kIM[rank_kIM == rank_cut_kIM]) %>%
      dplyr::mutate(is_in_final_module = rank_kIM <= rank_cut_kIM)
    
    pruning4 <- ggplot(pruned4, aes(x = rank_kIM, y = cumsum_kIM)) +
      geom_line() +
      geom_vline(aes(xintercept = rank_cut_kIM), linetype = "dashed", linewidth = 0.5) +
      geom_area(aes(fill = is_in_final_module), alpha = 0.3) +
      scale_fill_manual(values = c("#357D71", "#60AA98"), limits = c(T, F), guide = "none") +
      xlab(expression(paste("target rank based on ", italic(kIM)["pruned3"]))) +
      ylab(expression(paste("cumsum(", italic(kIM)["pruned3"], ")"))) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(size = 14, margin = margin(b = 0))) +
      labs(title = expression(paste(bold("Intermediate module after 3rd pruning "), "(", italic(M)["pruned3"], ")")),
           subtitle = paste0("Number of genes: ", nrow(pruned4)))
    
    if (n_iterations == 4) {
      
      final <- pruned4 %>%
        dplyr::filter(is_in_final_module)
      
      pruning5   <- ggplot() +
        labs(title = expression(paste(bold("Final module "), "(", italic(M)["final"], ")")),
             subtitle = paste0("Number of genes: ", nrow(final))) +
        theme_bw(base_size = 13) +
        theme(panel.background = element_blank(),
              plot.margin = margin(5.5, 46, -500, 5.5),
              plot.title = element_text(size = 14, margin = margin(b = 1)))
      
      return(pruning1  / plot_spacer() / pruning2  / plot_spacer() / pruning3  / plot_spacer() / pruning4 / plot_spacer() / pruning5 + plot_layout(heights = c(1, .08, 1, .08, 1, .08, 1, .08, 0.3)))
      
    } else {
      
      # pruning nr. 5
      pruned5 <- pruned4 %>%
        dplyr::filter(is_in_final_module) %>%
        arrange(desc(weight)) %>%
        dplyr::mutate(rank_IS = 1:length(target),
                      cumsum_IS = cumsum(weight),
                      rank_cut_IS = uik(rank_IS, cumsum_IS),
                      IS_cut = weight[rank_IS == rank_cut_IS]) %>%
        dplyr::mutate(is_in_final_module = rank_IS <= rank_cut_IS)
      
      pruning5 <- ggplot(pruned5, aes(x = rank_IS, y = cumsum_IS)) +
        geom_line() +
        geom_vline(aes(xintercept = rank_cut_IS), linetype = "dashed", linewidth = 0.5) +
        geom_area(aes(fill = is_in_final_module), alpha = 0.4) +
        scale_fill_manual(values = c("#06504A", "transparent"), limits = c(T, F), guide = "none") +
        new_scale_fill() +
        geom_area(aes(fill = weight > IS_cut, alpha = weight > IS_cut)) +
        scale_fill_manual(values = c("transparent", "#357D71"), limits = c(T, F), guide = "none") +
        scale_alpha_manual(values = c(0, 0.3), limits = c(T, F), guide = "none") +
        xlab(expression(paste("target rank based on ", italic(adj)["regulator"]))) +
        ylab(expression(paste("cumsum(", italic(adj)["regulator"], ")"))) +
        theme_bw(base_size = 13) +
        theme(plot.title = element_text(size = 14, margin = margin(b = 0.5))) +
        labs(title = expression(paste(bold("Intermediate module after 4th pruning "), "(", italic(M)["pruned4"], ")")),
             subtitle = paste0("Number of genes: ", nrow(pruned5)))
      
      final <- pruned5 %>%
        dplyr::filter(is_in_final_module)
      
      pruning6   <- ggplot() +
        labs(title = expression(paste(bold("Final module "), "(", italic(M)["final"], ")")),
             subtitle = paste0("Number of genes: ", nrow(final))) +
        theme_bw(base_size = 13) +
        theme(panel.background = element_blank(),
              plot.margin = margin(5.5, 46, -500, 5.5),
              plot.title = element_text(size = 14, margin = margin(b = 1)))
      
      return(pruning1  / plot_spacer() / pruning2  / plot_spacer() / pruning3  / plot_spacer() / pruning4 / plot_spacer() / pruning5 / plot_spacer() / pruning6 + plot_layout(heights = c(1, .08, 1, .08, 1, .08, 1, .08, 1, .08, 0.3)))
      
    }
    
  }
  
}


uncapitalize_first <- function(x) {
  
  first_char <- substr(x, 1, 1)
  second_char <- substr(x, 2, 2)
  
  if (grepl("[A-Z]", first_char) && !grepl("[A-Z]", second_char)) {
    # Only uncapitalize if second char is NOT uppercase
    return(paste0(tolower(first_char), substr(x, 2, nchar(x))))
  } else {
    return(x)
  }
}

wrapLongNames <- function(names, N, font_family = "Helvetica") {
  
  sapply(names, function(name) {
    
    grDevices::pdf(NULL)  # Opens a null PDF device for measuring text
    on.exit(grDevices::dev.off(), add = TRUE)
    
    name <- uncapitalize_first(name)
    
    visual_width <- strwidth(name, units = "inches", font = 1, family = font_family)
    
    if (visual_width <= N) {
      return(name)
    }
    
    chars <- unlist(strsplit(name, split = ""))
    partials <- sapply(seq_along(chars), function(i) {
      strwidth(paste(chars[1:i], collapse = ""), units = "inches", font = 1, family = font_family)
    })
    
    half_point_inch <- visual_width / 2
    
    # Find all space positions
    space_positions <- gregexpr(" ", name)[[1]]
    if (space_positions[1] == -1) {
      return(name)
    }
    
    # Find the space whose partial visual width is closest to the visual midpoint
    closest_space <- space_positions[which.min(abs(partials[space_positions] - half_point_inch))]
    
    # Insert linebreak at that space
    new_name <- paste0(
      substr(name, 1, closest_space - 1),
      "<br>",
      substr(name, closest_space + 1, nchar(name))
    )
    
    return(new_name)
    
  })
  
}





plotExprAlongPseudotime2 <- function(genes, sce, pseudotime_column = "pseudotime", cell_type_column = "cell_type", species_colors = NULL, cell_type_colors = NULL, font_size = 14, ncol = 1) {
  
  # if no species colors are provided, take the default
  species_names <- unique(sce$species)
  if (is.null(species_colors))
    species_colors <- species_color_ramp(seq(0, 1, len = length(species_names)))
  
  # extract metadata and logcounts of the genes of interest from the SCE object
  expr <- foreach::foreach(gene_name = genes,
                           .combine = dplyr::bind_rows,
                           .multicombine = TRUE) %do% {
                             
                             expr_gene <- data.frame(species = sce$species,
                                                     pseudotime = sce[[pseudotime_column]],
                                                     gene = gene_name,
                                                     expr = SingleCellExperiment::logcounts(sce)[gene_name,])
                             
                             if (!is.null(cell_type_column)) {
                               
                               expr_gene$cell_type <- sce[[cell_type_column]]
                               
                             }
                             
                             expr_gene
                             
                           } %>%
    dplyr::mutate(gene = factor(.data[["gene"]], levels = genes))
  
  # loess fit per gene and species
  loess_fit <- foreach::foreach(gene_name = genes,
                                .combine = dplyr::bind_rows,
                                .multicombine = TRUE) %:%
    foreach::foreach(species_name = species_names,
                     .combine = dplyr::bind_rows,
                     .multicombine = TRUE)  %do% {
                       
                       expr_gene_spec <- expr %>%
                         dplyr::filter(.data[["gene"]] == gene_name & .data[["species"]] == species_name)
                       fit  = stats::loess(expr ~ pseudotime, data = expr_gene_spec)
                       newx = data.frame(pseudotime = seq(min(expr_gene_spec$pseudotime), max(expr_gene_spec$pseudotime), len = 80))
                       pred = stats::predict(fit, newdata = newx, se = TRUE)
                       cbind(pseudotime = newx, expr = pred$fit, se = pred$se.fit) %>%
                         dplyr::mutate(moe = stats::qt(p = 0.975, df = nrow(expr_gene_spec) - 2)*.data[["se"]],
                                       lwr = .data[["expr"]] - .data[["moe"]],
                                       upr = .data[["expr"]] + .data[["moe"]],
                                       gene = gene_name,
                                       species = species_name)
                       
                     } %>%
    dplyr::mutate(gene = factor(.data[["gene"]], levels = genes))
  
  # initialize plot
  p <- expr %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["pseudotime"]], y = .data[["expr"]])) +
    ggplot2::xlab(pseudotime_column) +
    ggplot2::ylab("expression (logcounts)") +
    ggplot2::theme_bw(base_size = font_size)
  
  # if cell type information is present, add rug to the plot
  if (!is.null(cell_type_column)) {
    
    # if no cell type colors are provided, take the default
    if (is.null(cell_type_colors))
      cell_type_colors <- cell_type_color_ramp(seq(0, 1, len = length(unique(sce[[cell_type_column]]))))
    
    # calculate rug positions
    rug_positions <- loess_fit %>%
      dplyr::group_by(.data[["gene"]]) %>%
      dplyr::summarize(pos = min(.data[["lwr"]]) - 0.2*(max(.data[["upr"]]) - min(.data[["lwr"]]))) %>%
      tibble::deframe()
    
    # rug plot
    p <- p +
      ggplot2::geom_rug(data = expr %>% dplyr::mutate(expr = rug_positions[.data[["gene"]]]) %>% dplyr::filter(cell_type %in% c("Pluripotent_Cells", "Early_Ectoderm")),
                        ggplot2::aes(color = .data[["cell_type"]]),  sides = "b", show.legend = TRUE, length = ggplot2::unit(0.1, "npc"), linewidth = 0.15) +
      ggplot2::geom_rug(data = expr %>% dplyr::mutate(expr = rug_positions[.data[["gene"]]]) %>% dplyr::filter(cell_type == "Neurons"),
                        ggplot2::aes(color = .data[["cell_type"]]),  sides = "b", show.legend = TRUE, length = ggplot2::unit(0.1, "npc"), linewidth = 0.15) +
      ggplot2::scale_color_manual(values = cell_type_colors, labels = c("pluripotent\ncells", "early\nectoderm", "neurons"),  name = "cell type") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
      ggnewscale::new_scale_color() +
      theme(legend.key.spacing.y = unit(0.05, "cm"),
            legend.text = element_text(lineheight = 0.75, margin = margin(0.8, 5.5, 0.8, 5.5)))
    
  }
  
  # smooth curves
  p +
    ggplot2::geom_line(data = loess_fit, ggplot2::aes(color = .data[["species"]], group = .data[["species"]]), linewidth = 0.7) +
    ggplot2::geom_ribbon(data = loess_fit, ggplot2::aes(ymin = .data[["lwr"]], ymax = .data[["upr"]], fill = .data[["species"]]), alpha = 0.1) +
    ggplot2::scale_color_manual(values = species_colors) +
    ggplot2::scale_fill_manual(values = species_colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7), order = 1), fill = guide_legend(order = 1)) +
    ggplot2::facet_wrap(~.data[["gene"]], scales = "free_y", ncol = ncol, strip.position = "left") +
    ggplot2::scale_y_continuous(position = "right")
  
}


plotExprAlongPseudotimeAdjusted <- function(genes, sce, pseudotime_column = "pseudotime", cell_type_column = "cell_type", species_colors = NULL, cell_type_colors = NULL, font_size = 14, ncol = 1) {
  
  # if no species colors are provided, take the default
  species_names <- unique(sce$species)
  if (is.null(species_colors))
    species_colors <- species_color_ramp(seq(0, 1, len = length(species_names)))
  
  # extract metadata and logcounts of the genes of interest from the SCE object
  expr <- foreach::foreach(gene_name = genes,
                           .combine = dplyr::bind_rows,
                           .multicombine = TRUE) %do% {
                             
                             expr_gene <- data.frame(species = sce$species,
                                                     pseudotime = sce[[pseudotime_column]],
                                                     gene = gene_name,
                                                     expr = SingleCellExperiment::logcounts(sce)[gene_name,])
                             
                             if (!is.null(cell_type_column)) {
                               
                               expr_gene$cell_type <- sce[[cell_type_column]]
                               
                             }
                             
                             expr_gene
                             
                           } %>%
    dplyr::mutate(gene = factor(.data[["gene"]], levels = genes))
  
  # loess fit per gene and species
  loess_fit <- foreach::foreach(gene_name = genes,
                                .combine = dplyr::bind_rows,
                                .multicombine = TRUE) %:%
    foreach::foreach(species_name = species_names,
                     .combine = dplyr::bind_rows,
                     .multicombine = TRUE)  %do% {
                       
                       expr_gene_spec <- expr %>%
                         dplyr::filter(.data[["gene"]] == gene_name & .data[["species"]] == species_name)
                       fit  = stats::loess(expr ~ pseudotime, data = expr_gene_spec)
                       newx = data.frame(pseudotime = seq(min(expr_gene_spec$pseudotime), max(expr_gene_spec$pseudotime), len = 80))
                       pred = stats::predict(fit, newdata = newx, se = TRUE)
                       cbind(pseudotime = newx, expr = pred$fit, se = pred$se.fit) %>%
                         dplyr::mutate(moe = stats::qt(p = 0.975, df = nrow(expr_gene_spec) - 2)*.data[["se"]],
                                       lwr = .data[["expr"]] - .data[["moe"]],
                                       upr = .data[["expr"]] + .data[["moe"]],
                                       gene = gene_name,
                                       species = species_name)
                       
                     } %>%
    dplyr::mutate(gene = factor(.data[["gene"]], levels = genes))
  
  # initialize plot
  p <- expr %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["pseudotime"]], y = .data[["expr"]])) +
    ggplot2::xlab(pseudotime_column) +
    ggplot2::ylab("expression\n(logcounts)") +
    ggplot2::theme_bw(base_size = font_size)
  
  # if cell type information is present, add rug to the plot
  if (!is.null(cell_type_column)) {
    
    # if no cell type colors are provided, take the default
    if (is.null(cell_type_colors))
      cell_type_colors <- cell_type_color_ramp(seq(0, 1, len = length(unique(sce[[cell_type_column]]))))
    
    # calculate rug positions
    rug_positions <- loess_fit %>%
      dplyr::group_by(.data[["gene"]]) %>%
      dplyr::summarize(pos = min(.data[["lwr"]]) - 0.2*(max(.data[["upr"]]) - min(.data[["lwr"]]))) %>%
      tibble::deframe()
    
    # rug plot
    p <- p +
      ggplot2::geom_rug(data = expr %>% dplyr::mutate(expr = rug_positions[.data[["gene"]]]) %>% dplyr::filter(cell_type %in% c("Pluripotent_Cells", "Early_Ectoderm")),
                        ggplot2::aes(color = .data[["cell_type"]]),  sides = "b", show.legend = TRUE, length = ggplot2::unit(0.1, "npc"), linewidth = 0.15) +
      ggplot2::geom_rug(data = expr %>% dplyr::mutate(expr = rug_positions[.data[["gene"]]]) %>% dplyr::filter(cell_type == "Neurons"),
                        ggplot2::aes(color = .data[["cell_type"]]),  sides = "b", show.legend = TRUE, length = ggplot2::unit(0.1, "npc"), linewidth = 0.05) +
      ggplot2::scale_color_manual(values = cell_type_colors, labels = c("pluripotent\ncells", "early\nectoderm", "neurons"),  name = "cell type") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
      ggnewscale::new_scale_color() +
      theme(legend.key.spacing.y = unit(0.05, "cm"),
            legend.text = element_text(lineheight = 0.75, margin = margin(0.8, 5.5, 0.8, 5.5)),
            strip.text.y = element_text(size = 9.5),
            axis.text = element_text(size = 10))
    
  }
  
  # smooth curves
  p +
    ggplot2::geom_line(data = loess_fit, ggplot2::aes(color = .data[["species"]], group = .data[["species"]]), linewidth = 0.7) +
    ggplot2::geom_ribbon(data = loess_fit, ggplot2::aes(ymin = .data[["lwr"]], ymax = .data[["upr"]], fill = .data[["species"]]), alpha = 0.1) +
    ggplot2::scale_color_manual(values = species_colors) +
    ggplot2::scale_fill_manual(values = species_colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7), order = 1), fill = guide_legend(order = 1)) +
    ggplot2::facet_wrap(~.data[["gene"]], scales = "free_y", ncol = ncol, strip.position = "right")
  
}

species_color_ramp <- scales::colour_ramp(c("#07beb8", "steelblue3", "#2E4172", "#256F5C", "forestgreen", "#83AA3E", "#FFB600",  "#E28100"))

plotSmallTrees <- function(trees, species_colors = NULL, font_size = 16, tip_size = 1, branch_width = 0.4, ncol = NULL, out = "list") {
  
  if (inherits(trees, "phylo")) {
    
    if (is.null(species_colors))
      species_colors <- species_color_ramp(seq(0, 1, len = length(unique(trees$species))))
    
    tree_df <- getTreeDf(trees)
    
    return(
      plotSmallTree(trees, tree_df, species_colors, tip_size, branch_width, font_size)
    )
    
  } else {
    
    if (is.null(species_colors))
      species_colors <- species_color_ramp(seq(0, 1, len = length(unique(trees[[1]]$species))))
    
    tree_plots <- lapply(names(trees), function(name) {
      
      tree <- trees[[name]]
      tree_df <- getTreeDf(tree)
        
      plotSmallTree(tree, tree_df, species_colors, tip_size, branch_width, font_size)
      
    })
    
    names(tree_plots) <- names(trees)
    
    tree_plots <- tree_plots %>% matchScales(coord_flip = NULL, rev_x = NULL, rev_y = NULL)
    
    patchwork::wrap_plots(tree_plots, guides = "collect", ncol = ncol)
    
  }
  
}


plotSmallTree <- function(tree, tree_df, species_colors = NULL, tip_size = 1, branch_width = 0.4, font_size = 14) {
  
  suppressMessages(
      
      ggtree::ggtree(tree, layout = "daylight", size = branch_width) %<+%
        tree_df +
        ggtree::geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = tip_size*7, color = "transparent") +
        ggplot2::scale_fill_manual(values = species_colors, breaks = levels(tree_df$species)) +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = tip_size*5))) +
        ggplot2::theme(text = ggplot2::element_text(size = font_size),
                       legend.text = ggplot2::element_text(size = 0.8*22),
                       legend.title = ggplot2::element_text(size = 22),
                       legend.key.spacing.y = unit(0.15, "cm"),
                       plot.title = ggplot2::element_text(size = ggplot2::rel(1.2)))
      
    )
  
}

plotTreesAdjusted <- function(trees, species_colors = NULL, font_size = 14, tip_size = 1, branch_width = 0.4, rev_x = NULL, rev_y = NULL, coord_flip = NULL) {
  
  if (inherits(trees, "phylo")) {
    
    if (is.null(species_colors))
      species_colors <- species_color_ramp(seq(0, 1, len = length(unique(trees$species))))
    
    tree_df <- getTreeDf(trees)
    
    return(
      plotTree(trees, tree_df, species_colors, tip_size, branch_width, font_size)
    )
    
  } else {
    
    if (is.null(species_colors))
      species_colors <- species_color_ramp(seq(0, 1, len = length(unique(trees[[1]]$species))))
    
    tree_plots <- lapply(names(trees), function(name) {
      
      tree <- trees[[name]]
      tree_df <- getTreeDf(tree)
      
      if (name == "POU5F1") {
        
        plotTree_POU5F1(tree, tree_df, species_colors, tip_size, branch_width, font_size)
        
      } else if (name == "HOXA2") {
        
        plotTree_HOXA2(tree, tree_df, species_colors, tip_size, branch_width, font_size)
        
      } else if (name %in% coord_flip) {
        
        plotTree_coord_flip(tree, tree_df, species_colors, tip_size, branch_width, font_size)
        
      } else {
        
        plotTree(tree, tree_df, species_colors, tip_size, branch_width, font_size)
        
      }
      
    }) 
    names(tree_plots) <- names(trees)
    
    tree_plots %>% matchScales(rev_x = rev_x, rev_y = rev_y, coord_flip = coord_flip)
    
  }
  
}


plotTree <- function(tree, tree_df, species_colors = NULL, tip_size = 1, branch_width = 0.4, font_size = 14) {
    
    suppressMessages(
      
      ggtree::ggtree(tree, layout = "daylight", size = branch_width) %<+%
        tree_df +
        ggtree::geom_tiplab(ggplot2::aes(label = .data[["label"]], fill = .data[["species"]]), angle = 0, hjust = 0.5,  geom = "label", color = "black", size = tip_size, box.padding = ggplot2::unit(0.2, "lines"), label.padding = ggplot2::unit(0.2, "lines")) +
        ggplot2::scale_fill_manual(values = species_colors, breaks = levels(tree_df$species)) +
        ggplot2::theme(text = ggplot2::element_text(size = font_size),
                       legend.text = ggplot2::element_text(size = 0.8*14),
                       legend.title = ggplot2::element_text(size = 14, margin = margin(5.5, 5.5, 15, 0)),
                       legend.justification = "left",
                       legend.key.spacing.y = unit(0.05, "cm"),
                       legend.margin=margin(0,0,0,20),
                       plot.title = ggplot2::element_text(size = ggplot2::rel(1.2))) +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(color = "transparent", alpha = 1)))
      
    )
  
}

plotTree_coord_flip <- function(tree, tree_df, species_colors = NULL, tip_size = 1, branch_width = 0.4, font_size = 14) {
  
  suppressMessages(
    
    ggtree::ggtree(tree, layout = "daylight", size = branch_width) %<+%
      tree_df +
      ggtree::geom_tiplab(ggplot2::aes(label = .data[["label"]], fill = .data[["species"]]), angle = 0, hjust = 0.5,  geom = "label", color = "black", size = tip_size, box.padding = ggplot2::unit(0.2, "lines"), label.padding = ggplot2::unit(0.2, "lines")) +
      ggplot2::scale_fill_manual(values = species_colors, breaks = levels(tree_df$species)) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(color = "transparent"))) +
      ggplot2::theme(text = ggplot2::element_text(size = font_size),
                     legend.text = ggplot2::element_text(size = 0.8*22),
                     legend.title = ggplot2::element_text(size = 22),
                     legend.key.spacing.y = unit(0.15, "cm"),
                     plot.title = ggplot2::element_text(size = ggplot2::rel(1.2))) +
      coord_flip()
    
  )

}

plotTree_POU5F1 <- function(tree, tree_df, species_colors = NULL, tip_size = 1, branch_width = 0.4, font_size = 14) {
  
  suppressMessages(
      
      ggtree::ggtree(tree, layout = "daylight", size = branch_width) %<+%
        tree_df +
        geom_tiplab(data = . %>%
                      dplyr::filter(!(label %in% c("G1c1", "G1c2"))),
                    aes(label = label, fill = species), angle = 0, hjust = 0.5, geom = "label", color = "black", size = tip_size, box.padding = unit(0.2, "lines"), label.padding = unit(0.2, "lines")) +
        geom_tiplab(data = . %>%
                      dplyr::filter(label == "G1c1"),
                    aes(label = label, fill = species), angle = 0, hjust = 0, geom = "label", color = "black", size = tip_size, box.padding = unit(0.2, "lines"), label.padding = unit(0.2, "lines")) +
        geom_tiplab(data = . %>%
                      dplyr::filter(label == "G1c2"),
                    aes(label = label, fill = species), angle = 0, hjust = 0.5, vjust = 0, geom = "label", color = "black", size = tip_size, box.padding = unit(0.2, "lines"), label.padding = unit(0.2, "lines")) +
        ggplot2::scale_fill_manual(values = species_colors, breaks = levels(tree_df$species)) +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(color = "transparent"))) +
        ggplot2::theme(text = ggplot2::element_text(size = font_size),
                       legend.text = ggplot2::element_text(size = 0.8*14),
                       legend.title = ggplot2::element_text(size = 14),
                       legend.justification = "left",
                       legend.key.spacing.y = unit(0.15, "cm"),
                       plot.title = ggplot2::element_text(size = ggplot2::rel(1.2))) +
        coord_flip()
      
    )
  
}

plotTree_HOXA2 <- function(tree, tree_df, species_colors = NULL, tip_size = 1, branch_width = 0.4, font_size = 14) {
  
  suppressMessages(
      
      ggtree::ggtree(tree, layout = "daylight", size = branch_width) %<+%
        tree_df +
        geom_tiplab(data = . %>%
                      dplyr::filter(!(label %in% c("C2c2"))),
                    aes(label = label, fill = species), angle = 0, hjust = 0.5, geom = "label", color = "black", size = tip_size, box.padding = unit(0.2, "lines"), label.padding = unit(0.2, "lines")) +
        geom_tiplab(data = . %>%
                      dplyr::filter(label == "C2c2"),
                    aes(label = label, fill = species), angle = 0, hjust = 0.5, vjust = 0, geom = "label", color = "black", size = tip_size, box.padding = unit(0.2, "lines"), label.padding = unit(0.2, "lines")) +
        ggplot2::scale_fill_manual(values = species_colors, breaks = levels(tree_df$species)) +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(color = "transparent"))) +
        ggplot2::theme(text = ggplot2::element_text(size = font_size),
                       legend.text = ggplot2::element_text(size = 0.8*14),
                       legend.title = ggplot2::element_text(size = 14),
                       legend.justification = "left",
                       legend.key.spacing.y = unit(0.15, "cm"),
                       plot.title = ggplot2::element_text(size = ggplot2::rel(1.2)))
      
    )
  
}

matchScales <- function(plot_list, rev_x, rev_y, coord_flip) {
  
  x_ranges <- lapply(names(plot_list), function(name) {
    
    if (name %in% coord_flip) {
      
      c(ggplot2::layer_scales(plot_list[[name]])$y$range$range[1], ggplot2::layer_scales(plot_list[[name]])$y$range$range[2])
      
    } else {
      
      c(ggplot2::layer_scales(plot_list[[name]])$x$range$range[1], ggplot2::layer_scales(plot_list[[name]])$x$range$range[2])
      
    }
    
  })
  
  max_x_range <- max(sapply(x_ranges, function(vec) {
    
    vec[2] - vec[1]
    
  }))*1.1
  
  x_ranges_new <- lapply(x_ranges, function(vec) {
    
    x_range <- vec[2] - vec[1]
    x_expand <- (max_x_range - x_range) / 2
    c(vec[1] - x_expand, vec[2] + x_expand)
    
  })
  names(x_ranges_new) <- names(plot_list)
  
  y_ranges <- lapply(names(plot_list), function(name) {
    
    if (name %in% coord_flip) {
      
      c(ggplot2::layer_scales(plot_list[[name]])$x$range$range[1], ggplot2::layer_scales(plot_list[[name]])$x$range$range[2])
      
    } else {
      
      c(ggplot2::layer_scales(plot_list[[name]])$y$range$range[1], ggplot2::layer_scales(plot_list[[name]])$y$range$range[2])
      
    }
    
  })
  
  max_y_range <- max(sapply(y_ranges, function(vec) {
    
    vec[2] - vec[1]
    
  }))*1.1
  
  y_ranges_new <- lapply(y_ranges, function(vec) {
    
    y_range <- vec[2] - vec[1]
    y_expand <- (max_y_range - y_range) / 2
    c(vec[1] - y_expand, vec[2] + y_expand)
    
  })
  names(y_ranges_new) <- names(plot_list)
  
  plot_list <- lapply(names(plot_list), function(name) {
    
    if (name %in% rev_x & name %in% rev_y & name %in% coord_flip) {
      
      plot_list[[name]] +
        ggplot2::scale_x_continuous(limits = rev(y_ranges_new[[name]]), trans = "reverse") +
        ggplot2::scale_y_continuous(limits = rev(x_ranges_new[[name]]), trans = "reverse")
      
    } else if (name %in% rev_y & name %in% coord_flip) {
      
      plot_list[[name]] +
        ggplot2::scale_x_continuous(limits = rev(y_ranges_new[[name]]), trans = "reverse") +
        ggplot2::scale_y_continuous(limits = x_ranges_new[[name]])
      
    } else if (name %in% rev_x & name %in% coord_flip) {
      
      plot_list[[name]] +
        ggplot2::scale_x_continuous(limits = y_ranges_new[[name]]) +
        ggplot2::scale_y_continuous(limits = rev(x_ranges_new[[name]]), trans = "reverse")
      
    } else if (name %in% coord_flip) {
      
      plot_list[[name]] +
        ggplot2::scale_x_continuous(limits = y_ranges_new[[name]]) +
        ggplot2::scale_y_continuous(limits = x_ranges_new[[name]])
      
    } else if (name %in% rev_x & name %in% rev_y) {
      
      plot_list[[name]] +
        ggplot2::scale_x_continuous(limits = rev(x_ranges_new[[name]]), trans = "reverse") +
        ggplot2::scale_y_continuous(limits = rev(y_ranges_new[[name]]), trans = "reverse")
      
    } else if (name %in% rev_x) {
      
      plot_list[[name]] +
        ggplot2::scale_x_continuous(limits = rev(x_ranges_new[[name]]), trans = "reverse") +
        ggplot2::scale_y_continuous(limits = y_ranges_new[[name]])
      
    } else if (name %in% rev_y) {
      
      plot_list[[name]] +
        ggplot2::scale_x_continuous(limits = x_ranges_new[[name]]) +
        ggplot2::scale_y_continuous(limits = rev(y_ranges_new[[name]]), trans = "reverse")
      
    } else {
      
      plot_list[[name]] +
        ggplot2::scale_x_continuous(limits = x_ranges_new[[name]]) +
        ggplot2::scale_y_continuous(limits = y_ranges_new[[name]])
      
    }
    
  })
  
  plot_list
  
}


# helper function
plotDistMatsAdjusted <- function(dist_df, max_dist = 0.6, text_size = "small") {
  
  dist <- dist_df %>%
    dplyr::mutate(replicate1 = factor(replicate1, levels = unique(replicate1)),
                  replicate2 = factor(replicate2, levels = unique(replicate1)),
                  species1 = factor(species1, levels = unique(species1)),
                  species2 = factor(species2, levels = unique(species1)))
  
  dist <- dist %>%
    bind_rows(cbind(dist %>%
                      distinct(replicate1, species1),
                    dist %>%
                      distinct(replicate1, species1) %>%
                      dplyr::rename(replicate2 = replicate1, species2 = species1)))
  
  ggplot(dist, aes(x = replicate2, y = replicate1, fill = dist)) +
    geom_tile(colour="white", size=0.5) +
    scale_y_discrete(expand=c(0,0), limits = rev, position = "right") +
    scale_x_discrete(expand=c(0,0)) +
    scale_fill_gradientn(colors=rev(brewer.pal(name="RdYlGn", n=11)), limits= c(0, max_dist), breaks = seq(0, max_dist, by = 0.1),  name = "distance", na.value = "grey95") +
    theme_bw(base_size = ifelse(text_size == "small", 14, 17)) +
    xlab("replicate1") +
    ylab("replicate2") +
    theme(axis.text = element_text(size = case_when(text_size == "small" ~ 7.5, 
                                                    text_size == "medium" ~ 8.5,
                                                    text_size == "large" ~ 10.5)),
          axis.title = element_text(size = case_when(text_size == "small" ~ 14, 
                                                     text_size == "medium" ~ 16,
                                                     text_size == "large" ~ 16)),
          legend.margin=margin(0,0,0,20),
          legend.title = ggplot2::element_text(margin = margin(5.5, 5.5, 17, 0)),
          panel.spacing = unit(case_when(text_size == "small" ~ 0.1, 
                                         text_size == "medium" ~ 0.2,
                                         text_size == "large" ~ 0.2), "lines"),
          strip.text.x = element_text(margin = margin(0.12,0.12,0.12,0.12, "cm"), size = case_when(text_size == "small" ~ 10.8, 
                                                                                                   text_size == "medium" ~ 12.2,
                                                                                                   text_size == "large" ~ 17)),
          strip.text.y = element_text(margin = margin(0.12,0.12,0.12,0.12, "cm"), size = case_when(text_size == "small" ~ 10.8, 
                                                                                                   text_size == "medium" ~ 12.2,
                                                                                                   text_size == "large" ~ 17)),
          panel.border = element_rect(linewidth = 0.4),
          strip.background = element_rect(linewidth = 0.4),
          legend.justification = "left") +
    facet_grid(species1 ~ species2, scales = "free", space = "free", switch = "y")
  
}

# helper function to remove legend
no_legend <- function(p) {
  
  p + theme(legend.position = "none")
  
}

get_signif <- function(overlaps) {
  
  table <- overlaps %>% 
    pivot_wider(names_from = "in_module", values_from = "n") %>% 
    column_to_rownames("overlap_chipseq") %>% 
    as.matrix()
  
  fisher_test_result <- fisher.test(table)
  
  case_when(fisher_test_result$p.value < 0.001 ~ "***",
            fisher_test_result$p.value < 0.01 ~ "**",
            fisher_test_result$p.value < 0.05 ~ "*",
            T ~ "n.s.")
  
}


format_gencode_for_gviz <- function(gtf){
  
  seqlevelsStyle(gtf) <- "NCBI"
  
  gtf %>%
    as_tibble() %>%  
    filter(type %in% c("exon","UTR","CDS")) %>% 
    transmute(chromosome =seqnames, 
              start, end, width, strand, type,
              gene = gene_id,
              exon = exon_id,
              transcript = transcript_id,
              symbol = gene_name, gene_type) %>% 
    group_by(exon) %>% 
    dplyr::mutate( feature = case_when("CDS" %in% type ~ "protein_coding",
                                       "UTR" %in% type ~ "utr",
                                       T ~ gene_type)) %>% 
    dplyr::select( -gene_type, -type) %>% ungroup
  
}


format_gr_for_gviz <- function(gr){
  
  seqlevelsStyle(gr) <- "NCBI"
  
  gr %>%
    as_tibble() %>%
    dplyr::rename(chromosome = seqnames) %>% 
    dplyr::mutate(strand = "*")
  
}


plot_gviz <- function(gn,
                      gen,
                      gtf,
                      atac_bws,
                      atac_peaks,
                      extra_regions,
                      ltr7_hervh_gr,
                      tfbs_gr,
                      region,
                      chip_bws = NULL,
                      # chip_peaks = NULL,
                      max_atac = NA,
                      chrom_contact = NULL,
                      peak_labels) {
  
  # get chromosome
  chr <- as.character(as_tibble(region)$seqnames)
  start_region <- as_tibble(region)$start
  end_region <- as_tibble(region)$end
  
  # filter gene models for the region of interest
  gtf_filt <- gtf %>% 
    dplyr::filter(chromosome == chr &
                    ((start > start_region & start < end_region) | 
                       (end > start_region & end < end_region))) %>% 
    dplyr::rename(seqnames = chromosome) %>% 
    as_granges() %>% 
    group_by(gene, exon, transcript, symbol, feature) %>% 
    reduce_ranges_directed() %>% 
    as_tibble() %>% 
    dplyr::rename(chromosome = seqnames) %>% 
    dplyr::mutate(fill = ifelse(symbol == gn, "black", "grey60"),
                  font.color = ifelse(symbol == gn, "black", "grey60"))
  
  atac_peaks_filt <- atac_peaks %>% 
    dplyr::filter(chromosome == chr &
                    ((start > start_region & start < end_region) |
                       (end > start_region & end < end_region)))
  
  # if (!is.null(chip_peaks)) {
  #   
  #   chip_peaks_filt <- lapply(chip_peaks, function(x) {
  #     
  #     x %>% 
  #       dplyr::filter(chromosome == chr &
  #                       ((start > start_region & start < end_region) |
  #                          (end > start_region & end < end_region)))
  #     
  #   })
  #   
  # }
  
  ltr7_hervh_gr_filt <- ltr7_hervh_gr %>% 
    dplyr::filter(chromosome == chr &
                    ((start > start_region & start < end_region) |
                       (end > start_region & end < end_region)))
  
  ltr7_gr_filt <- ltr7_hervh_gr_filt %>% 
    dplyr::filter(type == "LTR7")
  
  # set correct chromosome format
  options(ucscChromosomeNames = F)
  
  # plot axis with genomic coordinates
  # idxTrack <- IdeogramTrack(genome=gen,  chromosome=chr)
  # idxTrack@bandTable$chrom <- gsub("chr", "", idxTrack@bandTable$chrom)
  # idxTrack@chromosome<-chr
  axis  <- GenomeAxisTrack(genome = gen)
  
  # plot GENCODE annotation
  gencode_track <- GeneRegionTrack(gtf_filt,
                                   chromosome = chr,
                                   genome = gen,
                                   showId = TRUE,
                                   geneSymbol = TRUE,
                                   col.group = gtf_filt$font.color,
                                   name= "GENCODE\nannotation",
                                   rotation.title=0,
                                   fontsize = 10, cex.title = 0.8,
                                   col.axis = "black", col.title = "black",
                                   fill = gtf_filt$fill,
                                   col = "transparent",
                                   collapseTranscripts = F, shape = "arrow")
  
  peak_colors <- c(hg38 = "#ACE293", gorGor6 = "#AAD0FF", macFas6 = "#e9b9dc")
  
  coverage_colors <- c(hg38 = "#4F904E", gorGor6 = "#416992", macFas6 = "#8804A8") 
  
  # plot ATAC-seq coverage 
  atac_coverage_tracks <- lapply(names(atac_bws), function(name) {
    
    DataTrack( range = atac_bws[[name]], 
               type = 'polygon', 
               chromosome = chr,
               name = paste0("ATAC-seq\n", name, " iPSCs"),
               fill.mountain = rep(coverage_colors[gen], 3),
               col.mountain = rep(coverage_colors[gen], 3),
               window = -1, windowSize = 100, genome = gen,
               cex.axis = 0.7,
               ylim = c(0, max_atac),
               col.title="black", cex.title=0.8,
               col.axis="black") 
    
  })
  
  # plot genome annotation info
  extra.track <- GeneRegionTrack(extra_regions,
                                 chromosome = chr,
                                 genome = gen,
                                 showId = TRUE,
                                 geneSymbol = TRUE,
                                 name= "genome\nalignment info",
                                 rotation.title=0,
                                 fontsize = 10, cex.title = 0.8,
                                 col.axis = "black", col.title = "black",
                                 fill = extra_regions$fill,
                                 col = "transparent",
                                 collapseTranscripts = F,
                                 stacking = "squish",
                                 just.group = "below",
                                 showOverplotting = TRUE)
  
  
  # plot LTR7 elements
  ltr7_track <- GeneRegionTrack(ltr7_hervh_gr_filt,
                                chromosome = chr,
                                genome = gen,
                                showId = TRUE,
                                geneSymbol = TRUE,
                                name= "LTR7 &\nHERVH-int",
                                rotation.title=0,
                                fontsize = 10, cex.title = 0.8,
                                fill = ltr7_hervh_gr_filt$fill, col = "transparent", col.axis = "black", col.title = "black", stacking = "pack",
                                just.group = "below",
                                showOverplotting = TRUE)
  
  # plot POU5F1 TFBS sites
  tfbs_track <- AnnotationTrack(tfbs_gr,
                                chromosome = chr,
                                genome = gen,
                                fill = tfbs_gr$fill,
                                name = "POU5F1 TFBS",
                                lwd = 0,
                                min.width = 1,
                                rotation.title=0,
                                fontsize = 10, cex.title = 0.8,
                                stacking = "dense",
                                showOverplotting = TRUE,
                                col = "transparent", col.axis = "black", col.title = "black")
  
  # plot chromatin interaction
  if (!is.null(chrom_contact)) {
    
    chrom_track <- GeneRegionTrack(chrom_contact,
                                   chromosome = chr,
                                   genome = gen,
                                   name= "contact w\npromoter",
                                   rotation.title=0,
                                   fontsize = 10, cex.title = 0.8,
                                   fill = "red4", col = "transparent", col.axis = "black", col.title = "black", stacking = "dense",
                                   showOverplotting = TRUE)
    
  }
  
  
  # plot LTR7 elements
  ltr7_track <- GeneRegionTrack(ltr7_hervh_gr_filt,
                                chromosome = chr,
                                genome = gen,
                                showId = TRUE,
                                geneSymbol = TRUE,
                                name= "LTR7 &\nHERVH-int",
                                rotation.title=0,
                                fontsize = 10, cex.title = 0.8,
                                fill = ltr7_hervh_gr_filt$fill, col = "transparent", col.axis = "black", col.title = "black", stacking = "pack",
                                just.group = "below",
                                showOverplotting = TRUE)
  
  # peak labels
  ortho.track <- GeneRegionTrack(peak_labels,
                                 chromosome = chr,
                                 genome = gen,
                                 showId = TRUE,
                                 geneSymbol = TRUE,
                                 rotation.title=0,
                                 name= "ortholgous CREs",
                                 fontsize = 10, cex.title = 0.8,
                                 col.axis = "black", col.title = "black",
                                 fill = "grey20",
                                 col = "transparent", 
                                 stacking = "squish",
                                 just.group = "below",
                                 showOverplotting = TRUE)
  
  
  chip_peak_colors <- c(hg38 = "darkseagreen3", macFas6 = "plum1")
  
  chip_coverage_colors <- c(hg38 = "#597759", macFas6 = "orchid4")
  
  # plot ChIP-seq coverage with peaks highlighted
  if (!is.null(chip_bws)) {
    
    chip.tracks.with.peaks <- lapply(names(chip_bws), function(name) {
      
      DataTrack( range = chip_bws[[name]],
                 type = 'polygon',
                 chromosome = chr,
                 name = paste0("POU5F1\nChIP-seq\n", name),
                 # rotation.title=0,
                 fill.mountain = rep(chip_coverage_colors[[gen]], 3),
                 col.mountain = rep(chip_coverage_colors[[gen]], 3),
                 window = -1, windowSize = 100, genome = gen,
                 col.title="black", cex.title=0.8,
                 col.axis="black",
                 ylim = c(0, 120),
                 yTicksAt = c(0, 40, 80, 120))
      
      # chip.peaks.track <- GeneRegionTrack(chip_peaks_filt[[name]],
      #                                     chromosome = chr,
      #                                     genome = gen,
      #                                     name = paste0("POU5F1\nChIP-seq\n", name),
      #                                     showId = F,
      #                                     geneSymbol = F,
      #                                     fill = chip_peak_colors[[gen]],
      #                                     col  = "transparent",
      #                                     col.title="black", cex.title=0.8,
      #                                     col.axis="black",
      #                                     shape = "box",
      #                                     stackHeight = 0.95)
      # 
      # 
      # OverlayTrack(trackList=list(chip.peaks.track, chip.track))
      
    })
    
  }
  
  track_list <- c(gencode_track, extra.track, ltr7_track, tfbs_track, atac_coverage_tracks)
  if (!is.null(chrom_contact)) track_list <- append(track_list, list(chrom_track), after = 3)
  if (!is.null(chip_bws)) track_list <- append(track_list, chip.tracks.with.peaks, after = length(track_list) - length(atac_coverage_tracks))
  
  # highlight LTR7 elements on all other tracks
  if (nrow(ltr7_gr_filt) + nrow(atac_peaks_filt) > 0) {
    
    tracks <- HighlightTrack(trackList = track_list,
                             start = c(atac_peaks_filt$start, ltr7_gr_filt$start),
                             end   = c(atac_peaks_filt$end, ltr7_gr_filt$end),
                             chromosome = chr,
                             alpha = 0.7,
                             fill = c(rep(peak_colors[gen], length(atac_peaks_filt$start)), rep("grey90", length(ltr7_gr_filt$start))),
                             col  = c(rep(peak_colors[gen], length(atac_peaks_filt$start)), rep("grey90", length(ltr7_gr_filt$start))))
    
    tracks <- list(tracks, ortho.track)
    
  } else {
    
    tracks <- list(track_list, ortho.track)
    
  }
  
  axis_size <- ifelse(gn == "SCGB3A2", 0.4, 0.45)
  
  extra_regions_reduced <- extra_regions %>% 
    dplyr::rename(seqnames = chromosome) %>% 
    as_granges() %>% 
    stretch(-2) %>% 
    reduce_ranges() %>% 
    as_tibble()
  
  extra_region_size <- case_when(nrow(extra_regions_reduced) < nrow(extra_regions) ~ 1.2,
                                 nrow(extra_regions) > 2 ~ 0.71,
                                 T ~ 0.59)
  
  gtf_size <- ifelse(length(unique(gtf_filt$transcript)) > 5, length(unique(gtf_filt$transcript))*0.15, 1.2)
  
  ortho_size = ifelse(gn == "SCGB3A2", 0.23, 0.3)
  
  # combine into a single plot
  track.list<- c(axis, tracks)
  plotTracks( track.list, 
              collapseTranscripts = F, shape = "arrow", 
              from = start_region, 
              to = end_region,
              title.width = 1.2,
              col.grid='grey' ,
              sizes = c(axis_size, gtf_size, extra_region_size, 0.47, rep(0.47, as.integer(!is.null(chrom_contact))), 0.47, rep(1, length(chip_bws)),  rep(1.1, length(atac_bws)), ortho_size),
              fontsize=11,
              cex.axis = 0.7,
              cex.main = 0.9,
              fontface.main = 1)
  
}