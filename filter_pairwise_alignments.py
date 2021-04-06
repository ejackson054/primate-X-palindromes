
"""
Filter regions of low-quality sequence alignment from pairwise DNA sequence alignments

Usage:  python filter_pairwise_alignments.py [infile] [outfile]

infile:  name of infile (alignment of two sequences in FASTA format)
outfile:  name of outfile (alignment of two sequences in FASTA format with low-quality regions masked with "N")

Default setting masks 100 bp windows with fewer than 60 matches; this setting can be changed by modifying "min_matches"

"""

from Bio import SeqIO
import numpy as np
import sys


"""Read in arguments"""
infile = sys.argv[1]
outfile = sys.argv[2]


"""Define functions"""

def alignment_to_dictionary(MSA):
    
    #Get bases for each position in the alignment:  "Dictionary of lists"
    MSA_len = len(MSA[0])
    bases_by_pos = {}
    
    for pos in np.arange(0,MSA_len):
        bases = []
        
        for sequence in MSA:
            letter = sequence.seq[pos]
            bases.append(letter[0])
                
        bases_by_pos[pos] = bases
      
    return(bases_by_pos)  
    
    
def get_fractions_by_window(step_start, window_size, bases_by_pos):
    
    not_applicable = 0
    same = 0
    different = 0
                
    for i in range(step_start, step_start + window_size):
                
        if ("-" in bases_by_pos[i]) or ("N" in bases_by_pos[i]):
            not_applicable += 1
    
        elif bases_by_pos[i][0] == bases_by_pos[i][1]:
            same += 1
            
        else:
            different+=1
        
    return(not_applicable, same, different)
    
    
    
"""Main loop"""
     
#Read in alignment
MSA = list(SeqIO.parse(infile, "fasta"))
MSA = [i.upper() for i in MSA]

#Get dictionary of bases by position
bases_by_pos = alignment_to_dictionary(MSA)
    
#Set minimum match (per step size)
min_matches = 60
        
#Set step size
window_size = 100
step_size = 1
number_steps = int(len(MSA[0].seq)/step_size) - window_size
step_start = 0
    
poor_MSA_starts = []
    
for step in range(1, number_steps + 1):

    na, same, different = get_fractions_by_window(step_start, window_size, bases_by_pos)
            
    if same < min_matches:
        poor_MSA_starts.append(step_start)
        
    step_start += step_size


#Convert list of poor window start sites (poor_MSA_starts) into list of all positions within those windows (poor_alignment_positions)
poor_alignment_positions = []
for j in poor_MSA_starts:
    positions = range(j, j+ window_size)
    for m in positions:
        if m not in poor_alignment_positions:
            poor_alignment_positions.append(m)
        
#Write out new alignment sequence
new_alignments = {}
for seq in MSA:
    new_alignments[seq.description] = ""
    for j in range(0,len(seq.seq)):
        if j not in poor_alignment_positions:
            new_alignments[seq.description] += seq.seq[j]
        else:
            new_alignments[seq.description] += "N"
            
#Write to new file
with open(outfile, "w") as new:
    for key in new_alignments.keys():
        new.write(">{0}\n{1}\n".format(key,new_alignments[key]))