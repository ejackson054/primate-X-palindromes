"""
Plots square dot plot annotated with positions of palindrome arms and protein-coding genes 

USAGE: python square_plots_annotated.py [square_plot] [amplicon_coordinates_s1] [gene_coordinates_s1] [sequence_s1] 
                                       [amplicon_coordinates_s2] [gene_coordinates_s2] [sequence_s2] [outfile]

    square_plot:                triangle dot plot created with fast_dot_plot.pl (see Sample Data, "CH251_397P16_671I19_prefinished_100_1_600.png")
    amplicon_coordinates:       tab-separated bedfile with amplicon coordinates (see Sample Data, "CH251_397P16_671I19_prefinished_amplicon_coordinates.bed")
    gene_coordinates:           tab-separated bedfile with gene coordinates (see Sample Data, "CH251_397P16_671I19_prefinished_genes.bed")
    sequence:                   FASTA sequence (see Sample Data, "CH251_397P16_671I19_prefinished.fasta")
    outfile:                    full path to save output pdf (e.g. "square_plot_annotated_P3_chimp_human.pdf")

"""

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import ndimage
import sys


"""Define variables"""

square_plot = sys.argv[1]
amp_coord_s1 = sys.argv[2]
gene_coord_s1 = sys.argv[3]
sequence_s1 = sys.argv[4]
amp_coord_s2 = sys.argv[5]
gene_coord_s2 = sys.argv[6]
sequence_s2 = sys.argv[7]
outfile = sys.argv[8]


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

              
def read_image_file(square_plot):

    im = plt.imread(square_plot)
    im = ndimage.rotate(im, 90)

    return(im)   
    

def extract_info(data, length1, length2, orientation):
    
    
    #Find which length is greater (so you can normalize arm position appropriately)
    if length1>length2:
        norm_factor = 500
        
    elif length2>=length1:
        norm_factor = 500 *length1/float(length2)
    
    
    #Find orientation (so you can normalize arm position appropriately)
    if orientation=="horizontal":
        
        def normalize(value):
            return(50 + (value/length1) * norm_factor)
        
    elif orientation == "vertical":
        
        def normalize(value):
            return(550 - (value/length1) * norm_factor)
    
    
    #Get raw values from table; normalize to length
    start_1 = normalize(float(data.loc[1,1]))
    end_1= normalize(float(data.loc[1,2]))
    start_2 = normalize(float(data.loc[3,1]))
    end_2 = normalize(float(data.loc[3,2]))
    length1_norm = normalize(length1)
        
    return(start_1,end_1, start_2, end_2, length1_norm)
        
   
                           
def plot_image(im,ax, length1, length2):

    ax.imshow(im, aspect="auto")
    
    ax.text(50, 520-(500*(length1/float(length2))),
                        "{0} x {1} kb".format(str(round(length1,-3))[0:-5],
                        str(round(length2,-3))[0:-5]),
                        fontsize=60, fontstyle="italic", fontweight="bold", ha="left")
        

def plot_on_axis_arms_horizontal(ax, start_1, end_1, start_2, end_2, length_1_norm):
    
    #Plot palindrome arrows
    ax.arrow(start_1, 560, end_1-start_1, 0, color='#1B75BC', linewidth=20, head_length = 3, head_width =7, length_includes_head = True, zorder=3)
    ax.arrow(start_2, 560, end_2-start_2, 0, color='#1B75BC', linewidth=20, head_length = 3, head_width = 7, length_includes_head = True, zorder=3)
            
    #Add thick #FBB040 lines
    ax.plot([start_1, start_2], [560,560], "-", color='#FBB040', linewidth=20, zorder = 1)
    ax.plot([50, end_1], [560,560], "-", color='black', linewidth=20, zorder = 1)
    ax.plot([end_2, length_1_norm], [560,560], "-", color='black', linewidth=20, zorder = 1)
        
               
def plot_on_axis_arms_vertical(ax, start_1, end_1, start_2, end_2, length_1_norm):
            
    #Plot palindrome arrows
    ax.arrow(40, start_1, 0, end_1-start_1, color='#1B75BC', linewidth=20, head_length = 3, head_width = 7, length_includes_head = True, zorder=3)
    ax.arrow(40, start_2, 0, end_2-start_2, color='#1B75BC', linewidth=20, head_length = 3, head_width = 7, length_includes_head = True, zorder=3)
            
    #Add thick #FBB040 lines
    ax.plot([40,40], [start_1, start_2], "-", color='#FBB040', linewidth=20, zorder = 1)
    ax.plot([40,40], [length_1_norm, end_2], "-", color='black', linewidth=20, zorder = 1)
    ax.plot([40,40], [end_1, 550], "-", color='black', linewidth=20, zorder = 1)
        
       
def plot_on_axis_genes_horizontal(ax, length1, length2, genes):
          
    if length1>length2:
        norm_factor = 500
        
    elif length2>=length1:
        norm_factor = 500 * length1/float(length2)
           
    #Plot gene arrows
    for i in genes.index:
    
        #Normalize gene lengths
        gene_start = 50 + (int(genes.loc[i,1])/length1) * norm_factor
        gene_end = 50 + (int(genes.loc[i,2])/length1) * norm_factor
        gene_name = genes.loc[i,4]

        ax.arrow(gene_start, 575, gene_end - gene_start,0, color='black', linewidth=3, head_length = 2, head_width = 6, length_includes_head = True)
        ax.text(gene_start, 590, gene_name, fontsize=fontsize, fontstyle="italic", rotation=60, ha="right", va="top") 

        
def plot_on_axis_genes_vertical(ax, length1, length2, genes):

    if length1>length2:
        norm_factor = 500
        
    elif length2>=length1:
        norm_factor = 500 * length1/float(length2)
        
    #Plot gene arrows
    for i in genes.index:
    
        #Normalize gene lengths
        gene_start = 550 - (int(genes.loc[i,1])/length1) * norm_factor
        gene_end = 550 - (int(genes.loc[i,2])/length1) * norm_factor
        gene_name = genes.loc[i,4]

        ax.arrow(30, gene_start, 0, gene_end - gene_start, color='black', linewidth=3, head_length = 2, head_width = 6, length_includes_head = True)
        ax.text(15, gene_start, gene_name, fontsize=fontsize, fontstyle="italic", rotation=60, ha="right", va="top") 
                    

"""Main code"""

#Set font size of labels on matplotlib plots
plt.rc('font', size=16)

#Set style of plots
sns.set_style('white')
fontsize = 60
 
#Read in arm coordinates
data_s1 = read_palindrome_file(amp_coord_s1)
data_s2 = read_palindrome_file(amp_coord_s2)

#Read in gene coordinates
genes_s1,length_s1 = read_gene_file(gene_coord_s1)
genes_s2,length_s2 = read_gene_file(gene_coord_s2)

#Read in images
im = read_image_file(square_plot)
   
#Get values for palindrome coordinates
#s1 = start1, s2 = start2, e1 = end1, e2 = end2
(s1_v, e1_v, s2_v, e2_v, length_norm_vert) = extract_info(data_s1,length_s1,length_s2,"vertical")
(s1_h, e1_h, s2_h, e2_h, length_norm_hor) = extract_info(data_s1,length_s1,length_s2,"horizontal")
     
#Initialize figure
fig = plt.figure(figsize = (50,50))

#Plot axes
ax = fig.add_axes([0, 1, 1, 1])

#Plot on Ax1
plot_image(im,ax, length_s1, length_s2)
       
#plot_on_axis_arms_vertical(ax, s1_v, e1_v, s2_v, e2_v, length_norm_vert)    
#plot_on_axis_arms_horizontal(ax, s1_h, e1_h, s2_h, e2_h, length_norm_hor)
    
plot_on_axis_genes_vertical(ax,length_s1,length_s2,genes_s1)    
plot_on_axis_genes_horizontal(ax,length_s1,length_s2, genes_s2)
    
#Remove axes
ax.set_axis_off()

fig.savefig(outfile,bbox_inches="tight",transparent=True)
    
        
