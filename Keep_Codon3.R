#set up & load packages
rm(list = ls())
options(stringsAsFactors = F)

suppressMessages(require(ape))
suppressMessages(require(getopt))


# Command Line Options
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

spec <- matrix(c(
  'help'       , 'h', 0, "logical"  , "Usage: Rscript Keep_Codon3.R input output",
  'input'      , 'f', 1, "character", "path to a .fasta nucleotide supermatrix",
  'output'     , 'o', 1, "character", "path to write the output .fasta nucleotide supermatrix with 3rd codon removed")
  , byrow = T, ncol = 5)

# Read options, give help/usage, error messages


opt <- getopt(spec)

if ( !is.null(opt$help) ){
  cat(getopt(spec, usage = T))
  q(status = 1)
}


if ( is.null(opt$input)     )  { stop("Input *.fasta  is required")                }
if ( is.null(opt$output)     )  { stop("Output path/name is required")                }

#

input_nt_fasta <- read.dna(opt$input, format="fasta")

Keep_Position_3 <- function(nucleotide_fasta){
  nucleotide_3R_fasta <- nucleotide_fasta[, seq(3, ncol(nucleotide_fasta), by = 3)]
  return(nucleotide_3R_fasta)
}

output_nt_codon3_only_fasta <- Keep_Position_3(input_nt_fasta)

write.FASTA(output_nt_codon3_only_fasta, opt$output)
message(paste("1st & 2nd Codon removed, and output written to",  opt$output))
  
