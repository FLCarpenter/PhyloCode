PhyloCode contains scripts designed for the preparation, processing, and manipulation of phylogenetic data, and the visualisation of phylogenetic trees. 
Some scripts are designed to be used in conjunction with https://github.com/tjcreedy/phylostuff. 
All scripts are designed for Command Line usage.


RY_Recoder.py Recodes nucleotide sequences from an input .fasta file replacing purines ("A", "a", "G", "g") with "R" and pyrimidines ("T", "t", "C", "c") with 'Y'. RY recoding aims to mitigate against the effect of base composition biases/heterogeneity in phylogenetic inference. 

Repartitioner.py is designed to expand on partitioner.py (see link above), it is designed specifically to be used with nexus formatted partitions produced by partitioner.py for nucleotide supermatrices. It takes an input partition and repartitions to match a .fasta with the 3rd codon position prior removed.
