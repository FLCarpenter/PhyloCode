PhyloCode contains scripts designed for the preparation, processing, and manipulation of phylogenetic data, and the visualisation of phylogenetic trees. 
Most scripts are designed to be used in conjunction with https://github.com/tjcreedy/pipelines. 
All scripts are designed for Command Line usage.


**RY_Recoder.py** DEPRECIATED - see below for updated recoder

**Nucleotide_Counter.py** DEPRECIATED - see below for updated counter



**Alignment_Counter.py** Counts each instance of nucleotides (A, T, G, C) or calculates the percentage of gaps (-) from an input .fasta file, giving an ouptut either per-sequence or alignment-wide results.
###### Usage
    python Nucleotide_Counter.py input.fasta [--gaps|-g] [--alignment|-a]
######
    python Alignment_Counter.py input.fasta
    python Alignment_Counter.py input.fasta -g
    python Alignment_Counter.py input.fasta -a
    python Alignment_Counter.py input.fasta -a -g
Per-sequence nucleotide counts:    *NO OPTIONAL ARGUMENTS*

Per-sequence gap percentage:      -g

Alignment-wide nucleotide counts: -a

Alignment-wide gap percentage:    -a -g


**RY_Recoder_Binary.py** Recodes nucleotide sequences from an input .fasta file replacing purines ("A", "a", "G", "g") with "0" or "R" and pyrimidines ("T", "t", "C", "c") with "1" or "Y". RY recoding aims to mitigate against the effect of base composition biases/heterogeneity in phylogenetic inference. 
###### Usage
    python RY_Recoder_Binary.py input.fasta input_datatype position output_datatype
###### 
    python RY_Recoder_Binary.py input.fasta nt N Binary
    python RY_Recoder_Binary.py input.fasta nt3r 1 RY

Note:
datatype: nt (nucleotide) or nt3r (nucleotide with 3rd codon position removed)

position: codons 1/2/3 or N (all)



**Remove_Codon3.R** is designed to allow for the removal of the 3rd codon since RAxML and other software don't support its removal through partitioning. The aligned input nucleotide supermatrix should in .fasta format and in the correct reading frame.
###### Usage
    Rscript Remove_Codon3.R --help|-h --input|-f input.fasta --output|-o output.fasta



**Repartitioner.py** is designed to expand on partitioner.py (see link above), it is designed specifically to be used with nexus formatted partitions produced by partitioner.py for nucleotide supermatrices. It takes an input partition and repartitions to match a .fasta with the 3rd codon position prior removed.
###### Usage
    python Repartitioner.py -f partition_type input.txt output.txt
###### 
    python Repartitioner.py -f gc input.txt output.txt
    python Repartitioner.py -f dc input.txt output.txt
    python Repartitioner.py -f c input.txt output.txt
(gc - g+codon12, dc - direction+codon12, c - codon12)
 



**raxml_partitioner.py** takes a nexus partition file and reformats it into a raxml compatible partition
###### Usage
    python raxml_partitioner.py input.nex output.txt



**split_supermatrix.py** splits a concatenated nucleotide alignment (supermatrix) into separate gene-specific FASTA files based on a user-provided partition file.
###### Usage
    python split_supermatrix.py input.txt input.fasta

input.txt: A partition file defining gene regions and their alignment coordinates.

input.fasta: A FASTA-formatted file containing the concatenated supermatrix alignment.


Note: Partition File Format - each line should specify a gene and its corresponding column range in the supermatrix, using the format:

path/to/gene_file.fasta = start-end


Example:

data/12S.fasta = 1-1699

data/16S.fasta = 1700-3200

