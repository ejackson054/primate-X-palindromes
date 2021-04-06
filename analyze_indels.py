"""
Find fraction of DNA sequence (Seq1) with no alignable counterpart in a second sequence (Seq2)

USAGE: python analyze_indels.py [bed_file] [amplicon_coordinates] [gene_coordinates] [sequence] [outfile]

    bed_file:              tab-separated bedfile created by bedgraph_from_alignments.py (see Sample Data, "Human_kmers_aligned_to_chimp_MAGED4.bed")
    amplicon_coordinates:  tab-separated bedfile with amplicon coordinates (see Sample Data, "CH251_397P16_671I19_prefinished_amplicon_coordinates.bed")
    gene_coordinates:      tab-separated bedfile with gene coordinates (see Sample Data, "CH251_397P16_671I19_prefinished_genes.bed")
    sequence:              FASTA sequence (see Sample Data, "CH251_397P16_671I19_prefinished.fasta")
    outfile:               name of text file containing results


"""

import pandas as pd
import itertools
import sys

"""Read in variables"""

bed_file = sys.argv[1]
amp_coord = sys.argv[2]
gene_coord = sys.argv[3]
sequence = sys.argv[4]
outfile = sys.argv[5]

"""Get fraction unaligned bases in each section"""
def get_fraction(unaligned_bases,all_bases):
    
    if len(all_bases)!=0:
        fraction = float(len(unaligned_bases))/len(all_bases)
    else:
        fraction = "None"
    
    return(fraction)
  
      
"""Get list of positions in problem areas"""
def get_problem_positions(sequence_file):
    
    #Get problem areas, if applicable
    with open(sequence_file,"r") as seq:
        for line in seq:
            if line[0]==">":
                try:
                    problem_areas = line.split("Problem areas: ")[1].split(",")
                    problem_areas = [i.strip() for i in problem_areas]
                    if problem_areas != "None":
                        problem_areas = [[int(i.split("-")[0]), int(i.split("-")[1])] for i in problem_areas]
                except:
                    problem_areas = ["None"]

    #GConvert list of problem areas into complete list of problem positions
    all_problem = []
    for i in problem_areas:
        if i!="None":
            problemrange = range(i[0],i[1])
            for j in problemrange:
                all_problem.append(j)
        
    return(all_problem)
    
    
"""Get fraction unaligned bases by section and gene status"""
def analyze_indels(bed_file, amplicon_coordinates_file, gene_coordinates_file, sequence_file):
    
    
    """Define indel positions"""
    
    #Get all "start" and "end" positions from bedfile as lists
    coverage_starts = []
    coverage_ends = []

    with open(bed_file,"r") as bed:
        for line in bed:
            coverage_starts.append(int(line.split("\t")[1]))
            coverage_ends.append(int(line.split("\t")[2]))
        
    #Find locations where >= 1000 bp do not align between species
    coverage_differences = []
    indel_coordinates = []

    for i in range(0,len(coverage_starts)-1):
        coverage_difference = coverage_starts[i+1] - coverage_ends[i] 
        if coverage_difference >= 1000:
            coverage_differences.append(coverage_difference)
            indel_coordinates.append([coverage_ends[i],coverage_starts[i+1]])
    
    """Read in gene coordinates"""
    
    gene_coordinates = pd.read_csv(gene_coordinates_file,sep="\t",header=None)
    
    #Get coordinates as lists of positions
    gene_coord = []
    geneless_coord = []
    
    for i in gene_coordinates.index:
        if gene_coordinates.ix[i,4]!="None":
            if gene_coordinates.ix[i,5]=="+":
                gene_coord.append([i for i in range(gene_coordinates.ix[i,1]-1000,gene_coordinates.ix[i,2])])
            elif gene_coordinates.ix[i,5]=="-":
                gene_coord.append([i for i in range(gene_coordinates.ix[i,1],gene_coordinates.ix[i,2]+1000)])
        else:
            geneless_coord.append([i for i in range(gene_coordinates.ix[i,1],gene_coordinates.ix[i,2])])
        
    #Combine these coordinates into a single range
    gene_coord = list(itertools.chain.from_iterable(gene_coord))
    geneless_coord = list(itertools.chain.from_iterable(geneless_coord))
       
       
    """Read in amplicon coordinates"""     
    #Read in amplicon coordinates
    amplicon_coordinates = pd.read_csv(amplicon_coordinates_file,sep="\t",header=None)

    #Define sections
    single_1 = list(range(int(amplicon_coordinates.ix[0,1])-1,int(amplicon_coordinates.ix[0,2])))
    single_2 = list(range(int(amplicon_coordinates.ix[4,1])-1,int(amplicon_coordinates.ix[4,2])))
    single = single_1 + single_2
    
    arm_1 = list(range(int(amplicon_coordinates.ix[1,1])-1,int(amplicon_coordinates.ix[1,2])))
    arm_2 = list(range(int(amplicon_coordinates.ix[3,1])-1,int(amplicon_coordinates.ix[3,2])))
    arms = arm_1 + arm_2
        
    spacer= list(range(int(amplicon_coordinates.ix[2,1])-1,int(amplicon_coordinates.ix[2,2])))


    """Find ALL bases by section"""
    
    #This returns a list of positions for every category
    single_genes_all = list(set(single).intersection(gene_coord))
    spacer_genes_all = list(set(spacer).intersection(gene_coord))
    arms_genes_all = list(set(arms).intersection(gene_coord))

    single_geneless_all = list(set(single).intersection(geneless_coord))
    spacer_geneless_all = list(set(spacer).intersection(geneless_coord))
    arms_geneless_all = list(set(arms).intersection(geneless_coord))
    

    """Find UNALIGNED bases by section"""
    all_unaligned = []
    
    for i in indel_coordinates:
        newrange = range(i[0],i[1])
        for j in newrange:
            all_unaligned.append(j)
    
    #This returns a list of positions for every category
    single_genes_unaligned = list(set(single).intersection(all_unaligned,gene_coord))
    spacer_genes_unaligned = list(set(spacer).intersection(all_unaligned,gene_coord))
    arms_genes_unaligned = list(set(arms).intersection(all_unaligned,gene_coord))
    
    single_geneless_unaligned = list(set(single).intersection(all_unaligned,geneless_coord))
    spacer_geneless_unaligned = list(set(spacer).intersection(all_unaligned,geneless_coord))
    arms_geneless_unaligned = list(set(arms).intersection(all_unaligned,geneless_coord))
    
    
    """Account for problem areas"""
    
    all_problem = get_problem_positions(sequence_file)
        
    single_genes_unaligned = list(set(single_genes_unaligned).difference(all_problem))
    spacer_genes_unaligned  = list(set(spacer_genes_unaligned).difference(all_problem))
    arms_genes_unaligned  = list(set(arms_genes_unaligned).difference(all_problem))
    
    single_geneless_unaligned  = list(set(single_geneless_unaligned).difference(all_problem))
    spacer_geneless_unaligned  = list(set(spacer_geneless_unaligned).difference(all_problem))
    arms_geneless_unaligned  = list(set(arms_geneless_unaligned).difference(all_problem))
        
    single_genes_all = list(set(single_genes_all).difference(all_problem))
    spacer_genes_all = list(set(spacer_genes_all).difference(all_problem))
    arms_genes_all = list(set(arms_genes_all).difference(all_problem))
        
    single_geneless_all= list(set(single_geneless_all).difference(all_problem))
    spacer_geneless_all = list(set(spacer_geneless_all).difference(all_problem))
    arms_geneless_all = list(set(arms_geneless_all).difference(all_problem))


    """Convert lists of positions into fractions"""
    fraction_single_genes_unaligned = get_fraction(single_genes_unaligned,single_genes_all)
    fraction_spacer_genes_unaligned = get_fraction(spacer_genes_unaligned,spacer_genes_all)
    fraction_arms_genes_unaligned = get_fraction(arms_genes_unaligned,arms_genes_all)

    fraction_single_geneless_unaligned =  get_fraction(single_geneless_unaligned,single_geneless_all)
    fraction_spacer_geneless_unaligned =  get_fraction(spacer_geneless_unaligned,spacer_geneless_all)
    fraction_arms_geneless_unaligned =  get_fraction(arms_geneless_unaligned,arms_geneless_all)
    
    """Return fraction unaligned for each category"""
    return(fraction_single_genes_unaligned, fraction_spacer_genes_unaligned,
    fraction_arms_genes_unaligned, fraction_single_geneless_unaligned,
    fraction_spacer_geneless_unaligned, fraction_arms_geneless_unaligned)
    
    
    
"""Main function"""

output = analyze_indels(bed_file, amp_coord, gene_coord, sequence)
   
with open(outfile, "w") as new:
    
    v1 = "Fraction unaligned (flanking, gene): {0}".format(output[0])
    v2 = "Fraction unaligned (spacer, gene): {0}".format(output[1])
    v3 = "Fraction unaligned (arm, gene): {0}".format(output[2])
    v4 = "Fraction unaligned (flanking, all other sequence): {0}".format(output[3])
    v5 = "Fraction unaligned (spacer, all other sequence): {0}".format(output[4])
    v6 = "Fraction unaligned (arm, all other sequence): {0}".format(output[5])
    
    new.write("{0}\n{1}\n{2}\n{3}\n{4}\n{5}".format(v1, v2, v3, v4, v5, v6))
             
  
