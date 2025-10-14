get_alignment_stats <- function(aligned_seqs, match = 1, mismatch = -1, gap_open = -5, gap_extend = -1) {
  if (length(aligned_seqs) != 2) {
    stop("Input must be a character vector of length 2: aligned reference and query sequences.")
  }
  
  ref <- aligned_seqs[1]
  qry <- aligned_seqs[2]
  
  if (nchar(ref) != nchar(qry)) {
    stop("Aligned sequences must be of the same length.")
  }
  
  ref_chars <- strsplit(ref, "")[[1]]
  qry_chars <- strsplit(qry, "")[[1]]
  
  len <- length(ref_chars)
  total_score <- 0
  matches <- 0
  mismatches <- 0
  aligned_pos <- 0
  gap_runs <- 0
  in_gap_ref <- FALSE
  in_gap_qry <- FALSE
  
  for (i in seq_along(ref_chars)) {
    a <- ref_chars[i]
    b <- qry_chars[i]
    
    # Check for gap in ref
    if (a == "-") {
      if (!in_gap_ref) {
        penalty <- if (i == 1 || i == len) gap_extend else gap_open
        total_score <- total_score + penalty
        gap_runs <- gap_runs + 1
        in_gap_ref <- TRUE
      } else {
        total_score <- total_score + gap_extend
      }
    } else {
      in_gap_ref <- FALSE
    }
    
    # Check for gap in qry
    if (b == "-") {
      if (!in_gap_qry) {
        penalty <- if (i == 1 || i == len) gap_extend else gap_open
        total_score <- total_score + penalty
        gap_runs <- gap_runs + 1
        in_gap_qry <- TRUE
      } else {
        total_score <- total_score + gap_extend
      }
    } else {
      in_gap_qry <- FALSE
    }
    
    # If either is a gap, skip scoring as substitution
    if (!a %in% AA_ALPHABET[1:20] | !b %in% AA_ALPHABET[1:20]) {
      next
    }
    
    aligned_pos <- aligned_pos + 1
      if (a == b) {
        matches <- matches + 1
        total_score <- total_score + match
      } else {
        mismatches <- mismatches + 1
        total_score <- total_score + mismatch
      }
  }
  
  xlen <- sum(ref_chars != "-")
  ylen <- sum(qry_chars != "-")
  al_prop <- aligned_pos / max(xlen, ylen)
  
  data.frame(
    xlen = xlen,
    ylen = ylen,
    al_length = aligned_pos,
    al_prop = al_prop,
    al_score = total_score,
    matches = matches,
    mismatches = mismatches,
    n_gaps = gap_runs,
    aa_cons = matches / (matches + mismatches)
  )
}