"""
Takes two FASTA sequences, and returns a bedgraph file describing the alignment of kmers from Sequence #2 onto Sequence #1 

USAGE: python Make_bedgraph_from_alignments_two_sequences.py [sequence1] [sequence2] [kmerlen] [alignment score] [outfile]

    sequence1:       Sequence file in FASTA format (sequence to which you are aligning)
    sequence2:       Sequence file in FASTA format (sequence to be broken into kmers)
    kmerLen:         100 for chimp versus human, 40 for rhesus versus human
    alignent score:  -12 for chimp versus human (up to 2 indels/mismatches), -18 for rhesus versus human  (up to 3 indels/mismatches)
    outfile:          e.g. Chimp_kmers_aligned_to_human_XAGE.bed

Purpose: Identifying indels between species
"""

from Bio import SeqIO
import os
import sys


"""Importing files and setting parameters"""

#Set arguments as variables
filename1 = sys.argv[1]
filename2 = sys.argv[2]
kmerLen = int(sys.argv[3])
alignScore = int(sys.argv[4])
outfile = sys.argv[5]

#Import sequence files
sequence1 = list(SeqIO.parse(filename1,"fasta"))[0].seq
sequence2 = list(SeqIO.parse(filename2,"fasta"))[0].seq

#Assign derivative names for future files  
kmer_file = "{0}_intermediate_file_with_kmers.fa".format(outfile.split(".")[0])
bowtie_index = outfile.split(".")[0]

samfile = "{}.sam".format(outfile.split(".")[0])
samfile_filtered = "{}_filtered.sam".format(outfile.split(".")[0])
bedgraph = outfile

#Set name of sequence being aligned to as "chrom"
with open(filename1,"r") as file1:
    for line in file1:
        if line[0]==">":
            chrom = line.split(">")[1].split(" ")[0].strip()


"""Getting kmers and aligning back to sequence file"""

#Make fasta file with all kmers from sequence2
with open(kmer_file,"w") as new:
    for i in range(len(sequence2)-kmerLen+1):
        new.write(">{0}\n{1}\n".format(i+1,sequence2[i:i+kmerLen]))
        
#Build index for sequence file
os.system("bowtie2-build {0} {1}".format(filename1,bowtie_index))

#Align all kmers back to sequence
os.system("n=8 bowtie2 -x {0} -f -U {1} -k 10 -S {2}".format(bowtie_index,kmer_file,samfile))

#Filter SAM file for only alignments with score >-11
with open(samfile,"r") as SAM:
    with open(samfile_filtered,"w") as new_sam:
        for line in SAM:
            if line[0]!="@":
                try:
                    score = int(line.split("AS:i:")[1].split("\t")[0])
                    if score >= alignScore:
                        new_sam.write("{0}".format(line))
                except:
                    pass

"""Counting depth at each position"""

#Count how many times each position in File 1 is covered by an alignment
with open(samfile_filtered,"r") as SAM:  
    
    #Initialize list of kmer positions
    positions = []
    
    #Read in every line and extract position
    for line in SAM:
        if line[0]!="@":
            position = int(line.split("\t")[3])
            positions.append(position)

#Sort results (slow step)
positions.sort()
counts = [positions.count(i) for i in positions]

#Group depths for bedgraph
initialize = 0
counter = 0
regions = []
depths = []

#Group all consecutive positions with the same depth
for i in range(len(positions)):
    if initialize==0:
        regions.append([positions[i]])
        depths.append([counts[i]])
        initialize+=1
    elif (counts[i]==depths[counter][-1]) and (positions[i]-positions[i-1]<2):
        regions[counter].append(positions[i])
        depths[counter].append(counts[i])
    else:
        regions.append([positions[i]])
        depths.append([counts[i]])
        counter+=1
        

#Write results to bedgraph
with open(bedgraph,"w") as bed:
    
    #Do this for every region
    for i in range(len(regions)):
        
        #Assign color based on depth
        if depths[i][0]==1:
            color = "0,0,0"       #black
        elif depths[i][0]==2:
            color = "0,255,0"     #green
        elif depths[i][0]>2:
            color = "255,0,0"     #red
        elif depths[i][0]==0:
            color = "255,255,255" #white
        
        #Assign region start
        if regions[i][0]==1:
            start = regions[i][0]
        else:
            start = regions[i][0]-1
        
        #Assign region end
        end = regions[i][-1]
        
        #Write to file
        bed.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(chrom,start,end,".",depths[i][0],".",start,end,color))

#Remove extra files at end
os.system("rm {}".format(kmer_file))
os.system("rm {}".format(samfile))
os.system("rm {}".format(samfile_filtered))
os.system("rm {}*.bt2".format(bowtie_index))
        
