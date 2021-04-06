"""Plots results from the file "fraction_alignable_pairwise.csv" (Fig. 4B)

USAGE: python fraction_alignable_pairwise_plot.py [infile] [outfile_pdf] [outfile_txt]

    infile:  full path to "fraction_alignable_pairwise.csv"
    outfile_pdf:  full path to save output pdf (e.g. "fraction_alignable_pairwise.pdf")
    outfile_txt:  full path to save output text results (e.g. "fraction_alignable_pairwise.txt")

infile should contain results from "fraction_alignable_pairwise.py" in the following format:
rownames = 12 conserved palindromes
colnames = results from all possible pairwise comparisons (single_human_vs_rhesus, arm_human_vs_chimp, spacer_human_vs_chimp, etc.)

"""

import scipy.stats as ss
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys


"""Set variables"""
infile = sys.argv[1]
outfile_pdf = sys.argv[2]
outfile_txt = sys.argv[3]


"""Reformat table for plotting"""
table = pd.read_table(infile, sep=",", index_col= "Unnamed: 0")
table.head()

values = []
types = []
species = []
full_types = []

for i in table.index:
    for j in table.columns:
        value = table.loc[i,j]
        if value != "None":
            values.append(float(value))
            
            if "chimp" in j.split("_"):
                species.append("chimp")
            elif "rhesus" in j.split("_"):
                species.append("rhesus")
            
            if "single" in j.split("_"):
                types.append("single")
            elif "arm" in j.split("_"):
                types.append("arm")
            elif "spacer" in j.split("_"):
                types.append("spacer")

df = pd.DataFrame()
df["Values"] = values
df["Types"] = types
df["Species"] = species

"""Plot results"""
fig, ax = plt.subplots()
palette = {"spacer": "#FBB040", "single":"silver","arm":"#1B75BC"}
fig = sns.barplot(data = df, x = "Species", y = "Values", hue = "Types", alpha = 0.6,
                  order = ["chimp", "rhesus"], palette = palette, saturation = 1)
fig = sns.swarmplot(data = df, x = "Species", y = "Values", hue = "Types",color = "black", size = 3,
                    order = ["chimp", "rhesus"], dodge = True)
fig.set_xticklabels(fig.get_xticklabels(), rotation=45, horizontalalignment='right')
fig.set_ylim(0,1.2)
ax.get_legend().remove()

plt.savefig(outfile_pdf, bbox_inches="tight",transparent=True)

"""Return averages and results of Mann-Whitney U"""

#Collect results:  Human versus chimp
single_human_vs_chimp = [float(i) for i in list(table["single_human_vs_chimp"]) if i !="None"]
single_chimp_vs_human = [float(i) for i in list(table["single_chimp_vs_human"]) if i !="None"]
single_chimp_human = single_human_vs_chimp + single_chimp_vs_human

spacer_human_vs_chimp = [float(i) for i in list(table["spacer_human_vs_chimp"]) if i !="None"]
spacer_chimp_vs_human = [float(i) for i in list(table["spacer_chimp_vs_human"]) if i !="None"]
spacer_chimp_human = spacer_human_vs_chimp + spacer_chimp_vs_human

arm_human_vs_chimp = [float(i) for i in list(table["arm_human_vs_chimp"]) if i !="None"]
arm_chimp_vs_human = [float(i) for i in list(table["arm_chimp_vs_human"]) if i !="None"]
arm_chimp_human = arm_human_vs_chimp  + arm_chimp_vs_human

#Collect results:  Human versus rhesus
single_human_vs_rhesus = [float(i) for i in list(table["single_human_vs_rhesus"]) if i !="None"]
single_rhesus_vs_human = [float(i) for i in list(table["single_rhesus_vs_human"]) if i !="None"]
single_rhesus_human = single_human_vs_rhesus + single_rhesus_vs_human

spacer_human_vs_rhesus = [float(i) for i in list(table["spacer_human_vs_rhesus"]) if i !="None"]
spacer_rhesus_vs_human = [float(i) for i in list(table["spacer_rhesus_vs_human"]) if i !="None"]
spacer_rhesus_human = spacer_human_vs_rhesus + spacer_rhesus_vs_human

arm_human_vs_rhesus = [float(i) for i in list(table["arm_human_vs_rhesus"]) if i !="None"]
arm_rhesus_vs_human = [float(i) for i in list(table["arm_rhesus_vs_human"]) if i !="None"]
arm_rhesus_human = arm_human_vs_rhesus  + arm_rhesus_vs_human

#Averages:  Human versus chimp
v1 = "Fraction alignable, human versus chimp (flanking): {0}".format(round(np.mean(single_chimp_human),3))
v2 = "Fraction alignable, human versus chimp (arms): {0}".format(round(np.mean(arm_chimp_human),3))
v3 = "Fraction alignable, human versus chimp (spacer): {0}".format(round(np.mean(spacer_chimp_human),3))

#Significance:  Human versus chimp
v4 = "Significance, human versus chimp (flanking vs. arms): {0}".format(ss.mannwhitneyu(arm_chimp_human,single_chimp_human))
v5 = "Significance, human versus chimp (flanking vs. spacer): {0}".format(ss.mannwhitneyu(spacer_chimp_human,single_chimp_human))
v6 = "Significance, human versus chimp (arms vs. spacer): {0}".format(ss.mannwhitneyu(spacer_chimp_human,arm_chimp_human))

#Averages:  Human versus rhesus
v7 = "Fraction alignable, human versus rhesus (flanking): {0}".format(round(np.mean(single_rhesus_human),3))
v8 = "Fraction alignable, human versus rhesus (arms): {0}".format(round(np.mean(arm_rhesus_human),3))
v9 = "Fraction alignable, human versus rhesus (spacer): {0}".format(round(np.mean(spacer_rhesus_human),3))

#Significance:  Human versus rhesus
v10 = "Significance, human versus chimp (flanking vs. arms): {0}".format(ss.mannwhitneyu(arm_rhesus_human,single_rhesus_human))
v11 = "Significance, human versus chimp (flanking vs. spacer): {0}".format(ss.mannwhitneyu(spacer_rhesus_human,single_rhesus_human))
v12 = "Significance, human versus chimp (arms vs spacer): {0}".format(ss.mannwhitneyu(spacer_rhesus_human,arm_rhesus_human))

#Write to file
with open(outfile_txt, "w") as new:
    new.write("{0}\n{1}\n{2}\n\n{3}\n{4}\n{5}\n\n{6}\n{7}\n{8}\n\n{9}\n{10}\n{11}".format(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12))



