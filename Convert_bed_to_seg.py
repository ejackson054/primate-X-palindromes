
"""Take a bedfile with depths for an entire chromosome, pull out region of interest, and convert to SEG file for visualization in IGV

USAGE: python Convert_bed_to_seg.py [IGV_sequence] [bedfile] [start]

IGV_sequence:  FASTA file to visualize in IGV
bedgraph:  bedGraph file with depths to visualize
start: Starting chromosome position for region of interest

"""


import pandas as pd
from Bio import SeqIO
import itertools
import sys


#Import arguments
igv = sys.argv[1]
bedgraph = sys.argv[2]
region_start = sys.argv[3]

#Import FASTA; get length, end position, name
seq = list(SeqIO.parse(igv,"fasta"))[0]
length = len(seq.seq)
region_end = region_start + length
name = seq.id

#Define functions
def import_bedgraph(bedgraph, name):
    
    chroms = []
    original_starts = []
    original_ends = []
    values = []
    
    with open(bedgraph, "r") as bed:

        for line in bed:
            fields = line.split("\t")
            start = int(fields[1])
            end = int(fields[2])
            depth = float(fields[3].split("\n")[0])
            if (start>=region_start) and (end<=region_end):
                original_starts.append(start)
                original_ends.append(end)
                chroms.append(name)
                if depth<0.5:
                    values.append(0)
                elif 0.5<=depth<=5:
                    values.append(2)
                elif depth>5:
                    values.append(4)

    return()

def find_empty_positions(df, length):
    
    #Get list of empty positions
    positions = range(1,length+1)
    ranges = [range(df.loc[i,"starts"],df.loc[i,"ends"]+1) for i in df.index]  
    ranges = list(set(itertools.chain.from_iterable(ranges)))
    not_in_bed =list(set(positions).difference(ranges)) 
    
    #Sort into groups
    groups = []
    count = 0
    for i in not_in_bed:
        if groups == []:
            groups.append([i])
        elif i==(groups[count][-1]+1):
            groups[count].append(i)
        else:
            groups.append([i])
            count+=1
    
    return(groups)
    
    
def main(bedgraph,igv,start):
    
    outfile = "{0}.seg".format(igv.split(".fasta")[0])

    #Import bedgraph
    chroms, original_starts, original_ends, values = import_bedgraph(bedgraph)

    #Subtract down starts and ends
    new_starts = [i-original_starts[0]+1 for i in original_starts]
    new_ends =  [i-original_starts[0]+1 for i in original_ends]
    
    #Combine into table
    df = pd.DataFrame(columns=["sequence_names","chroms","starts","ends","values"])
    df["sequence_names"] =  [name for i in new_starts]
    df["starts"] = new_starts
    df["ends"] = new_ends
    df["chroms"] = chroms
    df["values"] = values
    
    #Find coordinates NOT in bed
    groups = find_empty_positions(df, length)

    #Make these into new entries
    sequence_names_2 = []
    starts_2 = []
    ends_2 = []
    chroms_2 = []
    values_2 = []

    for i in groups:
        sequence_names_2.append(name)
        chroms_2.append(name)
        starts_2.append(i[0])
        ends_2.append(i[-1])
        values_2.append(0)
    
    #Make this into table
    df2 = pd.DataFrame(columns=["sequence_names","chroms","starts","ends","values"])
    df2["sequence_names"]=sequence_names_2
    df2["chroms"]=chroms_2
    df2["starts"]=starts_2
    df2["ends"]=ends_2
    df2["values"]=values_2
    
    #Merge tables
    df3 = pd.concat([df,df2])    
    df3 = df3.sort_values(["starts"])
    df3["new_index"]=range(0,len(df3.index))
    df3 = df3.set_index("new_index")
    
    #Write to output
    with open(outfile,"w") as new:
        for i in df3.index:
            sn = df3.loc[i,"sequence_names"]
            chrom = df3.loc[i,"chroms"]
            start = df3.loc[i,"starts"]
            end = df3.loc[i,"ends"]
            value = df3.loc[i,"values"]
            new.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(sn,chrom,start,end,value))
    

if __name__ == '__main__':
    main()
