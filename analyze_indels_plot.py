"""
Find fraction of DNA sequence (Seq1) with no alignable counterpart in a second sequence (Seq2)

USAGE: python analyze_indels_plot.py [base_dir] [outfile_pdf] [outfile_txt]

    base_dir:      directory containing six .csv files containing results from analyze_indels.py (e.g. see Sample Data, "analyze_indels_human_aligned_to_rhesus.csv")
    outfile_pdf:   path to save output pdf (e.g. "analyze_indels.pdf")
    outfile_txt:   path to save output txt (e.g. "analyze_indels.txt")

"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import scipy.stats as ss
import numpy as np
import sys


"""Define variables"""
base_dir = sys.argv[1]
outfile_pdf = sys.argv[2]
outfile_txt = sys.argv[3]

os.chdir(base_dir)


"""Define function"""
#Source:  https://kite.com/python/examples/702/scipy-compute-a-confidence-interval-from-a-dataset

def get_95_ci(values):
    
    n = len(values)
    std_err = ss.sem(values)
    h = std_err * ss.t.ppf((1 + 0.95) / 2, n - 1)
    return(h)


"""Read in tables"""

h_to_r = pd.read_csv("Analyze_indels_three_way_conserved_genes_human_aligned_to_rhesus.csv",index_col="Unnamed: 0")
r_to_h = pd.read_csv("Analyze_indels_three_way_conserved_genes_rhesus_aligned_to_human.csv",index_col="Unnamed: 0")
h_to_c = pd.read_csv("Analyze_indels_three_way_conserved_genes_human_aligned_to_chimp.csv",index_col="Unnamed: 0")
c_to_h = pd.read_csv("Analyze_indels_three_way_conserved_genes_chimp_aligned_to_human.csv",index_col="Unnamed: 0")
c_to_r = pd.read_csv("Analyze_indels_three_way_conserved_genes_chimp_aligned_to_rhesus.csv",index_col="Unnamed: 0")
r_to_c = pd.read_csv("Analyze_indels_three_way_conserved_genes_rhesus_aligned_to_chimp.csv",index_col="Unnamed: 0")


"""Reformat tables into single long-format dataframe"""

#Find total length of values
tables_list = [h_to_r, r_to_h ,h_to_c, c_to_h, c_to_r, r_to_c]
species_list = ["Human_aligned_to_rhesus","Rhesus_aligned_to_human","Human_aligned_to_chimp","Chimp_aligned_to_human",
           "Chimp_aligned_to_rhesus","Rhesus_aligned_to_chimp"]
total_len = (len(h_to_r.index) + len(r_to_h.index) + len(h_to_c.index) + len(c_to_h.index) + len(c_to_r.index) + len(r_to_c.index)) * 6

#Modify label names
dictionary = {"single":"Flanking","spacer":"Spacer","arms":"Arms"}

#Initialize dataframe
tidy_df = pd.DataFrame(index=range(0,total_len),columns=["Sequence type","Palindrome", "Sequence region","Species","Fraction unaligned"])

values = []
seq_types = []
seq_regions = []
palindromes = []
species_combos = []

for number in range(0,len(tables_list)):
    
    table = tables_list[number]
    species = species_list[number]

    for j in table.index:
        for k in table.columns:
            palindromes.append(j)
            values.append(table.loc[j,k])
            seq_types.append(str(k).split("_")[-1])
            seq_regions.append(str(k).split("_")[-2])
            species_combos.append(species)

tidy_df["Species"] = species_combos
tidy_df["Palindrome"] = palindromes
tidy_df["Sequence type"] = seq_types
tidy_df["Sequence region"] = seq_regions
tidy_df["Fraction unaligned"] = values
tidy_df = tidy_df[tidy_df["Fraction unaligned"]!="None"]
tidy_df["Fraction unaligned"] = [float(i) for i in tidy_df["Fraction unaligned"]]


"""Plot figure"""
fig = plt.figure()
ax = fig.add_subplot(111)

x = [1,2,3.5,4.5,6,7]

single_gene = tidy_df[tidy_df["Sequence region"]=="single"]
single_gene = single_gene[single_gene["Sequence type"]=="gene"]
single_gene_mean = np.mean(list(single_gene["Fraction unaligned"]))
single_gene_std = get_95_ci(list(single_gene["Fraction unaligned"]))

single_geneless = tidy_df[tidy_df["Sequence region"]=="single"]
single_geneless = single_geneless[single_geneless["Sequence type"]=="geneless"]
single_geneless_mean = np.mean(list(single_geneless["Fraction unaligned"]))
single_geneless_std = get_95_ci(list(single_geneless["Fraction unaligned"]))

arm_gene = tidy_df[tidy_df["Sequence region"]=="arms"]
arm_gene = arm_gene[arm_gene["Sequence type"]=="gene"]
arm_gene_mean = np.mean(list(arm_gene["Fraction unaligned"]))
arm_gene_std = get_95_ci(list(arm_gene["Fraction unaligned"]))

arm_geneless = tidy_df[tidy_df["Sequence region"]=="arms"]
arm_geneless = arm_geneless[arm_geneless["Sequence type"]=="geneless"]
arm_geneless_mean = np.mean(list(arm_geneless["Fraction unaligned"]))
arm_geneless_std = get_95_ci(list(arm_geneless["Fraction unaligned"]))

spacer_gene = tidy_df[tidy_df["Sequence region"]=="spacer"]
spacer_gene = spacer_gene[spacer_gene["Sequence type"]=="gene"]
spacer_gene_mean = np.mean(list(spacer_gene["Fraction unaligned"]))
spacer_gene_std = get_95_ci(list(spacer_gene["Fraction unaligned"]))

spacer_geneless = tidy_df[tidy_df["Sequence region"]=="spacer"]
spacer_geneless = spacer_geneless[spacer_geneless["Sequence type"]=="geneless"]
spacer_geneless_mean = np.mean(list(spacer_geneless["Fraction unaligned"]))
spacer_geneless_std = get_95_ci(list(spacer_geneless["Fraction unaligned"]))

y = [single_gene_mean, single_geneless_mean, arm_gene_mean, arm_geneless_mean, spacer_gene_mean, spacer_geneless_mean]

ax.bar(x[0], y[0], width = 1, color = "black", yerr = single_gene_std)
ax.bar(x[1], y[1], width = 1, color = "black", alpha = 0.6, yerr = single_geneless_std, hatch = "//")
ax.bar(x[2], y[2], width = 1, color = "#1B75BC", yerr = arm_gene_std)
ax.bar(x[3], y[3], width = 1, color = "#1B75BC", alpha = 0.6, yerr = arm_geneless_std,hatch = "//")
ax.bar(x[4], y[4], width = 1, color = "#FBB040", yerr = spacer_geneless_std)
ax.bar(x[5], y[5], width = 1, color = "#FBB040", alpha = 0.6, yerr = spacer_geneless_std,hatch = "//")

ax.set_ylim(0,0.35)
plt.savefig(outfile_pdf, format="pdf")


"""Run statistics"""

#Subset dataframe by type:  single_gene, single_geneless, etc.
single = tidy_df[tidy_df["Sequence region"]=="single"]
single_gene = single[single["Sequence type"]=="gene"]
single_geneless = single[single["Sequence type"]=="geneless"]

arm = tidy_df[tidy_df["Sequence region"]=="arms"]
arm_gene = arm[arm["Sequence type"]=="gene"]
arm_geneless = arm[arm["Sequence type"]=="geneless"]

spacer = tidy_df[tidy_df["Sequence region"]=="spacer"]
spacer_gene = spacer[spacer["Sequence type"]=="gene"]
spacer_geneless = spacer[spacer["Sequence type"]=="geneless"]

#Calculate means
spacer_gene_mean = np.mean(list(spacer_gene["Fraction unaligned"]))
spacer_geneless_mean = np.mean(list(spacer_geneless["Fraction unaligned"]))
arm_gene_mean = np.mean(list(arm_gene["Fraction unaligned"]))
arm_geneless_mean = np.mean(list(arm_geneless["Fraction unaligned"]))
single_gene_mean = np.mean(list(single_gene["Fraction unaligned"]))
single_geneless_mean = np.mean(list(single_geneless["Fraction unaligned"]))

#Calculate statistics
single_sig = ss.mannwhitneyu(single_gene["Fraction unaligned"], single_geneless["Fraction unaligned"])
arm_sig = ss.mannwhitneyu(arm_gene["Fraction unaligned"], arm_geneless["Fraction unaligned"])
spacer_sig = ss.mannwhitneyu(spacer_gene["Fraction unaligned"], spacer_geneless["Fraction unaligned"])

#Write to file
with open(outfile_txt, "w") as new:
    
    v1 = "Average (spacer, gene): {0}".format(round(spacer_gene_mean,4))
    v2 = "Average (spacer, all other sequence): {0}".format(round(spacer_geneless_mean,4))
    v3 = "Average (arm, gene): {0}".format(round(arm_gene_mean,4))
    v4 = "Average (arm, all other sequence): {0}".format(round(arm_geneless_mean,4))
    v5 = "Average (flanking, gene): {0}".format(round(single_gene_mean,4))
    v6 = "Average (flanking, all other sequence): {0}".format(round(single_geneless_mean,4))
    v7 = "Significance (spacer, gene versus all other sequence): {0}".format(spacer_sig)
    v8 = "Significance (arm, gene versus all other sequence): {0}".format(arm_sig)
    v9 = "Significance (flanking, gene versus all other sequence): {0}".format(single_sig)

    new.write("{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n\n{6}\n{7}\n{8}".format(v1,v2,v3,v4,v5,v6,v7,v8,v9))