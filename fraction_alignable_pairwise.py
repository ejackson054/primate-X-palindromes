"""
Takes a pairwise alignment and returns the fraction sequence that can be aligned relative to each sequence

USAGE: python fraction_alignable_pairwise.py [sequence_1] [sequence_2]
    
    sequence_1: raw pairwise alignment generated with ClustalW (.fasta format)
    sequence_2: filtered pairwise alignment after processing with filter_pairwise_alignments.py (.fasta format)

"""


from Bio import SeqIO
import sys

unfiltered_alignment = sys.argv[1]
filtered_alignment = sys.argv[2]

def get_fraction_alignable(unfiltered_alignment, filtered_alignment):
                
    unfilt = list(SeqIO.parse(unfiltered_alignment, "fasta"))
    filt = list(SeqIO.parse(filtered_alignment, "fasta"))
            
    seq1_len_unfilt = len([i for i in unfilt[0].seq if i not in ["-","N"]])
    seq1_len_filt = len([i for i in filt[0].seq if i not in ["-","N"]])

    seq2_len_unfilt = len([i for i in unfilt[1].seq if i not in ["-","N"]])
    seq2_len_filt = len([i for i in filt[1].seq if i not in ["-","N"]])
            
    fraction_seq1_alignable = float(seq1_len_filt)/seq1_len_unfilt
    fraction_seq2_alignable = float(seq2_len_filt)/seq2_len_unfilt
            
    return(fraction_seq1_alignable, fraction_seq2_alignable)
            
    
fraction_seq1_alignable, fraction_seq2_alignable  = get_fraction_alignable(unfiltered_alignment, filtered_alignment)

print("Fraction seq1 alignable: {0}".format(fraction_seq1_alignable))
print("Fraction seq1 alignable: {0}".format(fraction_seq2_alignable))