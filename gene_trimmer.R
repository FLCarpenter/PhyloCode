#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(getopt)
})


# ----------------------------
# Define command-line options
# ----------------------------
spec <- matrix(c(
  'input',      'i', 1, "character", "Codon-aligned nucleotide FASTA (single mt PCG)",
  'output',     'o', 1, "character", "Output FASTA file or directory",
  'threshold',  't', 2, "double",    "AA occupancy threshold [default = 0.5]"
), ncol = 5, byrow = TRUE)

opts <- getopt(spec)

# ----------------------------
# Check required options
# ----------------------------
if (is.null(opts$input) || is.null(opts$output)) {
  cat(getopt(spec, usage=TRUE))
  stop("ERROR: --input and --output must be specified")
}

threshold <- ifelse(is.null(opts$threshold), 0.5, opts$threshold)

# --------------------------------------------------
# Read aligned nucleotide FASTA
# --------------------------------------------------
read_fasta <- function(fasta_file) {
  read.dna(fasta_file, format = "fasta")
}

# --------------------------------------------------
# Translate nucleotide alignment to amino acids
# --------------------------------------------------
translate_alignment <- function(dna, genetic_code = 5) {
  trans(dna, code = genetic_code)
}

# --------------------------------------------------
# Find AA trimming bounds based on occupancy threshold
# --------------------------------------------------
find_trim_bounds <- function(aa, threshold = 0.5) {
  
  # Convert AAbin to character matrix
  aa_char <- matrix(NA, nrow = nrow(aa), ncol = ncol(aa))
  for (i in seq_len(nrow(aa))) {
    aa_char[i, ] <- as.character(aa[i, ])
  }
  
  # Occupancy: proportion of real residues (not "-" or "X")
  occupancy <- colMeans(aa_char != "-" & aa_char != "X")
  
  valid_sites <- which(occupancy > threshold)
  
  if (length(valid_sites) == 0) {
    stop("ERROR: No amino-acid sites exceed the occupancy threshold.")
  }
  
  left <- valid_sites[1]
  right <- tail(valid_sites, 1)
  
  if (left >= right) {
    stop("ERROR: Invalid trimming bounds detected.")
  }
  
  list(left = left, right = right)
}


# --------------------------------------------------
# Trim nucleotide alignment using AA bounds
# --------------------------------------------------
trim_nt_alignment <- function(dna, bounds) {
  
  nt_start <- (bounds$left - 1) * 3 + 1
  nt_end <- bounds$right * 3
  
  dna[, nt_start:nt_end]
}

# --------------------------------------------------
# Remove sequences with internal stop codons
# --------------------------------------------------
remove_stop_codons <- function(dna, genetic_code = 5) {
  
  aa <- trans(dna, numcode = genetic_code)
  aa_mat <- as.matrix(aa)
  
  keep <- apply(aa_mat, 1, function(x) !any(x == "*"))
  
  dna[keep, ]
}

# --------------------------------------------------
# Resolve output path and create directories if needed
# --------------------------------------------------
resolve_output_path <- function(output, input) {
  
  # If output is a directory or ends with "/"
  if (grepl("/$", output) || dir.exists(output)) {
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    return(file.path(output, basename(input)))
  }
  
  # Otherwise output is a file path
  out_dir <- dirname(output)
  if (out_dir != ".") {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  output
}

# ----------------------------
# Main workflow
# ----------------------------
dna <- read_fasta(opts$input)
aa <- translate_alignment(dna)

bounds <- find_trim_bounds(aa, threshold)
dna_trimmed <- trim_nt_alignment(dna, bounds)

out_path <- resolve_output_path(opts$output, opts$input)
  
write.FASTA(dna_trimmed, file = out_path)

cat("Trimmed alignment written to:", out_path, "\n")
