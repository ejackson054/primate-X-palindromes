"""Identify tissue-specific genes from a log2-normalized expression matrix using the tissue specificity index (TSI)


Usage: python Tissue_specific_genes.py [expression_table] [specificity]

expression_table: CSV file containing gene expression by tissue (transcripts per million, TPM)
specificity: Minimum specificity for declaring gene "tissue-specific" (TSI must be between 0 and 1)


For overview of TSI, see "A benchmark of gene expression tissue-specificity metrics," Kryuchkova-Mostacci and Robinson-Rechavi 2016

"""

import numpy as np
import pandas as pd
import sys

#Set variables
expression_table = sys.argv[1]  
anno_file = sys.argv[2] 
specificity_threshold = sys.argv[3] 

#Import expression table
dataframe = pd.read_table(expression_table, sep=",", index_col="Genes")

#Remove genes expressed <1 TPM in all tissues
cutoff = np.log2(2)
dataframe = dataframe[(dataframe>cutoff).any(axis=1)]

#Calculate TSI
dataframe_sum = np.sum(dataframe,axis=1)
dataframe_specificity = dataframe.div(dataframe_sum,axis=0)

#Identify tissue-specific genes 
genes = []
tissues = []
values = []

for i in dataframe_specificity.index:
    for j in dataframe_specificity.columns:         
        if (dataframe_specificity.loc[i,j]>=specificity_threshold) and (dataframe.loc[i,j] == np.max(dataframe.loc[i,:])):
            genes.append(i)
            tissues.append(j)
            values.append(dataframe_specificity.loc[i,j])

#Write to CSV file with gene, chromosome, tissue of specificity, TSI
genes = pd.Series(genes)
values = pd.Series(values)
tissues = pd.Series(tissues)
dictionary = {"Genes":genes,"Tissues":tissues,"Values":values}
df = pd.DataFrame(dictionary)
df = df.set_index("Genes")
df.to_csv("Tissue_specific_genes.csv")


