
sample_species <- function(target_df, bin_counts, n_total, seed = 100) {
  set.seed(seed)
  bins <- names(bin_counts)
  target_df$used <- FALSE  # track which rows have been sampled
  sampled_list <- list()
  
  # 1st pass: sample per bin_counts as much as possible
  for (bin in bins) {
    in_bin <- target_df[target_df$bin_pt == bin & !target_df$used, ]
    desired_n <- bin_counts[[bin]]
    available_n <- nrow(in_bin)
    n <- min(desired_n, available_n)
    
    if (n > 0) {
      sampled <- in_bin[sample(seq_len(n), n), , drop = FALSE]
      target_df[rownames(sampled), "used"] <- TRUE
      sampled_list[[bin]] <- sampled
    } else {
      sampled_list[[bin]] <- NULL
    }
  }
  
  # Combine all initial samples
  all_sampled <- do.call(rbind, sampled_list) 
  all_sampled <- all_sampled[!is.na(all_sampled$pseudotime),]
  
  # If we still need more samples to reach n_total, sample from the rest
  remaining_needed <- n_total - nrow(all_sampled)
  
  if (remaining_needed > 0) {
    remaining_pool <- target_df[!target_df$used, ]
    if (nrow(remaining_pool) >= remaining_needed) {
      fill <- remaining_pool[sample(seq_len(nrow(remaining_pool)), remaining_needed), , drop = FALSE]
    } else {
      warning("Not enough cells in target_df to reach desired sample size.")
      fill <- remaining_pool  # take all we can
    }
    all_sampled <- rbind(all_sampled, fill)
  }
  
  # Clean up
  all_sampled$used <- NULL
  return(all_sampled)
}


run_nonrandom_sampling<-function(coldata, bin, n_quantile_bins = 15, seed){
  
  cyno<-coldata[coldata$bin==bin & coldata$species=="cynomolgus",]
  hum<-coldata[coldata$bin==bin & coldata$species=="human",]
  gor<-coldata[coldata$bin==bin & coldata$species=="gorilla",]
  
  # identify the minimum one
  datasets <- list(cyno = cyno, hum = hum, gor = gor)
  smallest_name <- names(which.min(sapply(datasets, nrow)))
  # Access the dataset with the fewest rows
  smallest_dataset <- datasets[[smallest_name]]
  
  print(paste("smallest dataset is from", smallest_dataset$species[1], 
              "with", nrow(smallest_dataset), "rows"))
  
  # Quantile-based breaks from the minimal dataset
  quantile_breaks <- quantile(smallest_dataset$pseudotime, 
                              probs = seq(0, 1, length.out = n_quantile_bins + 1), 
                              na.rm = TRUE)
  
  quantile_breaks <- quantile_breaks[!duplicated(quantile_breaks)]
  print(paste("number of effective breaks", length(quantile_breaks)))
  
  # Cut all species using cyno-based quantile bins
  smallest_dataset$bin_pt  <- cut(smallest_dataset$pseudotime,  
                                  breaks = quantile_breaks, include.lowest = TRUE)
  # Count how many samples per bin in the smallest dataset
  bin_counts <- table(smallest_dataset$bin_pt)
  
  # count also in the other tables
  cyno$bin_pt <- cut(cyno$pseudotime, breaks = quantile_breaks, include.lowest = TRUE)
  hum$bin_pt  <- cut(hum$pseudotime,  breaks = quantile_breaks, include.lowest = TRUE)
  gor$bin_pt  <- cut(gor$pseudotime,  breaks = quantile_breaks, include.lowest = TRUE)
  
  # Cyno sampling
  print("sampling cyno")
  sample_cyno <- sample_species(
    target_df = cyno %>% rownames_to_column("BC"),
    bin_counts = bin_counts,
    n_total = nrow(smallest_dataset),
    seed = seed
  )
  
  # Human sampling
  print("sampling human")
  sample_hum <- sample_species(
    target_df = hum %>% rownames_to_column("BC"),
    bin_counts = bin_counts,
    n_total = nrow(smallest_dataset),
    seed = seed
  )
  
  # Gorilla sampling
  print("sampling gorilla")
  sample_gor <- sample_species(
    target_df = gor %>% rownames_to_column("BC"),
    bin_counts = bin_counts,
    n_total = nrow(smallest_dataset),
    seed = seed
  )
  return(bind_rows(sample_cyno, sample_hum, sample_gor))
}


