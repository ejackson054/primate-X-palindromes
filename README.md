# primate-X-palindromes



Supplemental Code for Jackson et al. 2021, "Natural selection preserves large palindromes on the primate X Chromosome"

Written by H. Skaletsky (fast_dot_plot.pl, Perl) and E. Jackson (all other codes, Python 3)

For scripts used for analysis of human palindrome deletions from the 1000 Genomes Project, see https://github.com/lsteitz/y-amplicon-evolution (Teitz et al. 2018 ). 

"Sample Data" contains example files used as input for these scripts (see description within each script for reference to appropriate sample data)



**fast_dot_plot.pl** :  Creates square and triangular dot plots from FASTA files

**triangle_plot_annotated.py**:  Annotates triangle dot plot with positions of palindrome arms and protein-coding genes

**square_plot_annotated.py**:  Annotates square dot plot with positions of palindrome arms and protein-coding genes

**filter_pairwise_alignments.py**:  Masks low-quality regions from pairwise FASTA sequence alignments

**tissue_specific_genes.py**:  Finds tissue-specific genes from a CSV containing gene expression levels (TPM) averaged by tissue

**fraction_alignable_pairwise.py**:  Finds fraction of each sequence from a pairwise FASTA sequence alignment with high-quality alignment to the other sequence.  General metric for sequence orthology.

**fraction_alignable_pairwise_plot.py**:  Plots results collected from fraction_alignable_pairwise.py across all palindromes and species comparisons

**analyze_indels.py**:  Finds fraction of each sequence from a pairwise FASTA sequence without alignable counterpart (unaligned stretches of >1kb) in the other sequence. Detects large insertions/deletions.

**analyze_indels_plot.py**:  Plot results collected from analyze_indels.py aross all palindromes and species comparisons

**bedgraph_from_alignments.py**:  Takes two FASTA sequences, and returns a bedgraph file describing the alignment of kmers from one sequence onto the other

**convert_bed_to_seg.py**:  Takes bedfile with depths for an entire chromosome, pulls out region of interest, and converts to SEG file for visualization in IGV
