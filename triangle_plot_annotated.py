# -*- coding: utf-8 -*-
"""
Plots triangle dot plot annotated with positions of palindrome arms and protein-coding genes 

USAGE: python triangle_plots_annotated.py [triangle_plot] [amplicon_coordinates] [gene_coordinates] [sequence] [outfile]

    triangle_plot:              triangle dot plot created with fast_dot_plot.pl (see Sample Data, "CH251_397P16_671I19_prefinished_100_1_600.png")
    amplicon_coordinates:       tab-separated bedfile with amplicon coordinates (see Sample Data, "CH251_397P16_671I19_prefinished_amplicon_coordinates.bed")
    gene_coordinates:           tab-separated bedfile with gene coordinates (see Sample Data, "CH251_397P16_671I19_prefinished_genes.bed")
    sequence:                   FASTA sequence (see Sample Data, "CH251_397P16_671I19_prefinished.fasta")
    outfile:                    full path to save output pdf (e.g. "triangle_plot_annotated_P3_chimp.pdf")

"""

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys


"""Define variables"""
triangle_plot = sys.argv[1]
amp_coord = sys.argv[2]
gene_coord = sys.argv[3]
sequence = sys.argv[4]
outfile = sys.argv[5]


"""Define functions"""

def read_palindrome_file(amp_coord):
    
    data = pd.read_csv(amp_coord,sep="\t",header=None)
    arm_1_start = data.loc[1,1]
    arm_1_end = data.loc[1,2]
    
    data.loc[1,1] = arm_1_end
    data.loc[1,2] = arm_1_start
   
    return(data)
        

def read_gene_file(gene_coord):
    
    genes = pd.read_csv(gene_coord, sep="\t", header=None)
    length =  float(genes.loc[(max(genes.index)),2])  
    genes = genes.loc[[i for i in genes.index if genes.loc[i,4]!="None"],:]
        
    #For genes that are in reverse orientation:  Flip the start and end positions
    for i in genes.index:
        if genes.loc[i,3]=="-":
            start = genes.loc[i,1]
            end = genes.loc[i,2]
            genes.loc[i,1] = end
            genes.loc[i,2] = start
            genes.loc[i,6] = end
            genes.loc[i,7] = start    
          
    return(genes, length)

        
def read_image_file(triangle_plot):
    
    im = plt.imread(triangle_plot)   
    return(im)   
    

def extract_info(data):
        
    #Get raw values
    start_1_raw = float(data.loc[1,1])
    end_1_raw = float(data.loc[1,2])
    start_2_raw = float(data.loc[3,1])
    end_2_raw = float(data.loc[3,2])
    length = float(data.loc[4,2])
        
    #Normalize to length
    start_1 = 50 + (1-(length - start_1_raw)/length) * 500
    end_1= 50 + (1-(length - end_1_raw)/length) * 500
    start_2 = 50 + (1-(length - start_2_raw)/length) * 500
    end_2 = 50 + (1-(length - end_2_raw)/length) * 500     
        
    return(start_1,end_1,start_2,end_2)
        

def plot_on_axis_arms(ax, length, start_1, end_1, start_2, end_2, im):
    
    #Plot image
    ax.imshow(im,aspect="auto")
              
    #Plot palindrome arrows
    ax.arrow(start_1, 560,end_1-start_1,0, color='#1B75BC', linewidth=8,head_length = 3, head_width = 7, length_includes_head = True, zorder=3)
    ax.arrow(start_2, 560,end_2-start_2, 0, color='#1B75BC', linewidth=8,head_length = 3, head_width = 7, length_includes_head = True, zorder=3)
            
    #Add lines for spacer and flanking
    ax.plot([start_1, start_2],[560,560],"-",color='#FBB040', linewidth=8, zorder = 2)
    ax.plot([50, end_1],[560,560],"-",color='black', linewidth=8, zorder = 1)
    ax.plot([end_2, 550],[560,560],"-",color='black', linewidth=8, zorder = 1)
                         
def plot_on_axis_genes(ax, length, genes):

    #Plot gene arrows
    for i in genes.index:
    
        #Normalize gene lengths
        gene_start = 50 + (1-(length - int(genes.loc[i,1]))/length) * 500
        gene_end = 50 + (1-(length - int(genes.loc[i,2]))/length) * 500
        gene_name = genes.loc[i,4].strip()

        #Plot arrows and names
        ax.arrow(gene_start, 575,gene_end - gene_start,0, color='black', linewidth=3,head_length = 2, head_width = 6, length_includes_head = True)
        ax.text(gene_start, 590, gene_name, fontsize=24,fontstyle="italic",rotation=60,ha="right", va="top") 
                                                                                                                                                                                                                     

"""Main code"""

#Set font size of labels on matplotlib plots
plt.rc('font', size=16)

#Set style of plots
sns.set_style('white')

#Read in arm coordinates
data = read_palindrome_file(amp_coord)

#Read in gene coordinates 
genes, length = read_gene_file(gene_coord)
    
#Read in images
im = read_image_file(triangle_plot)

#Get values for palindrome coordinates
(start_1, end_1, start_2, end_2) = extract_info(data)
    
#Initialize figure
fig = plt.figure(figsize = (20,20))
    
#Plot axes
ax = fig.add_axes([0, 300, 1,1])

#Plot on Axes
plot_on_axis_arms(ax, length, start_1, end_1, start_2, end_2, im)
plot_on_axis_genes(ax, length, genes)
    
#Remove axes
ax.set_axis_off()

fig.savefig(outfile, bbox_inches="tight",transparent=True)