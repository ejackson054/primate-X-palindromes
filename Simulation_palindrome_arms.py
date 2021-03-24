"""Simulate evolution of ancestral palindrome with ongoing gene conversion between arms with variable degrees of GC bias"""

import numpy as np
from scipy.stats import binom
import sys
import os



"""Set parameters"""

iteration = int(sys.argv[1])
GC_content = float(sys.argv[2])   #Median palindrome arms: 0.46 (from Supplemental Table 1 GC)
GC_bias = float(sys.argv[3])

length = 37000                    #median human palindrome length, rounded to nearest kb:  for all palindromes, 37000
fraction_differences = 0.00047    #median of all 36 conserved X palindromes in human, chimp and rhesus

relative_frequencies = {"AT_to_TA": 0.115, "AT_to_CG": 0.136, "AT_to_GC": 0.566,
                        "CG_to_AT": 0.227, "CG_to_TA":0.833, "CG_to_GC": 0.210}

substitution_rate_branch34 = 0.0000000120
conversion_rate_branch34 = (2*substitution_rate_branch34) / fraction_differences     #both of these rates are events per nucleotide, per generation

substitution_rate_branch1 = 0.0000000186
conversion_rate_branch1 = (2*substitution_rate_branch1) / fraction_differences     #both of these rates are events per nucleotide, per generation

substitution_rate_branch2 = 0.0000000207
conversion_rate_branch2 = (2*substitution_rate_branch2) / fraction_differences     #both of these rates are events per nucleotide, per generation

AT_content = 1 - GC_content
AT_bias = 1 - GC_bias

output_dir = "/lab/solexa_page/emily/GC_content_simulation/02_04_2021/simulations/simulated_arms/GC_bias_{0}_GC_content_{1}".format(GC_bias, GC_content)

if not os.path.isdir(output_dir):
    os.system("mkdir {0}".format(output_dir))
    
    

"""Define functions"""

def separate_sequence_by_base(sequence):
    
    positions_A = [i for i in range(0,len(sequence)) if sequence[i]=="A"]
    positions_T = [i for i in range(0,len(sequence)) if sequence[i]=="T"]
    positions_C = [i for i in range(0,len(sequence)) if sequence[i]=="C"]
    positions_G = [i for i in range(0,len(sequence)) if sequence[i]=="G"]

    position_dict = {"A": positions_A, "T": positions_T, "C": positions_C, "G": positions_G}
    
    return(position_dict)
    

def calculate_substitution_rates_AT_CG(relative_frequencies, substitution_rate):
    
    #Absolute frequencies
    absolute_frequencies = {}
    
    for i in relative_frequencies.keys():
        absolute_frequencies[i] = relative_frequencies[i]*substitution_rate
    
    #AT and GC rates
    AT_rate = absolute_frequencies["AT_to_TA"] + absolute_frequencies["AT_to_GC"] + absolute_frequencies["AT_to_CG"]
    CG_rate = absolute_frequencies["CG_to_GC"] + absolute_frequencies["CG_to_AT"] + absolute_frequencies["CG_to_TA"]
    
    return(AT_rate, CG_rate)
    
def generate_substitution_matrix(relative_frequencies):
    
    sum_AT = relative_frequencies["AT_to_TA"] + relative_frequencies["AT_to_GC"] + relative_frequencies["AT_to_CG"]
    sum_CG = relative_frequencies["CG_to_GC"] + relative_frequencies["CG_to_AT"] + relative_frequencies["CG_to_TA"]
    
    #Specific rates
    AT_to_TA = relative_frequencies["AT_to_TA"]/sum_AT
    AT_to_GC = relative_frequencies["AT_to_GC"]/sum_AT
    AT_to_CG = relative_frequencies["AT_to_CG"]/sum_AT
    CG_to_GC = relative_frequencies["CG_to_GC"]/sum_CG
    CG_to_AT = relative_frequencies["CG_to_AT"]/sum_CG
    CG_to_TA = relative_frequencies["CG_to_TA"]/sum_CG
    
    #Convert to matrix
    substitution_matrix = {"A": [0, AT_to_TA, AT_to_CG, AT_to_GC], 
                "T": [AT_to_TA, 0, AT_to_CG, AT_to_GC], 
                "G": [CG_to_AT, CG_to_TA, 0, CG_to_GC], 
                "C": [CG_to_AT, CG_to_TA, CG_to_GC, 0]}
    
    return(substitution_matrix)


    
def assign_conversion_prob(base_1, base_2):
    
    c1 = (base_1 in ["G","C"]) and (base_2 in ["G","C"])
    c2 = (base_1 in ["A","T"]) and (base_2 in ["A","T"])
    c3 = (base_1 in ["G","C"]) and  (base_2 in ["A","T"])
    c4  = (base_1 in ["A","T"]) and (base_2 in ["G","C"])
    
    if c1 or c2:
        
        prob = [0.5, 0.5]
        
    elif c3:
        
        prob = [GC_bias, AT_bias]
    
    elif c4:
        
        prob = [AT_bias, GC_bias]
    
    return(prob)
        
    

def mutate_arm(arm, possible_positions, number_mut, position_dict):
    
    global substitution_matrix
    
    pos = np.random.choice(possible_positions, size=number_mut)
    
    for i in pos:
        old = arm[i]
        new = np.random.choice(bases,p = substitution_matrix[old])
        arm[i] = new
        position_dict[old].remove(i)
        position_dict[new].append(i)
        
  
def convert_arm(arm_1, arm_2, number_conv, position_dict_1, position_dict_2):  
    
    possible_pos = list(range(0,len(arm_1)))
    
    pos = np.random.choice(possible_pos, size=number_conv)
        
    for i in pos:
        
        prob = assign_conversion_prob(arm_1[i], arm_2[i])
        
        converted_base = np.random.choice([arm_1[i], arm_2[i]],p = prob)
        
        #Update list of positions that contains each base!
        if converted_base == arm_1[i]:
            position_dict_2[arm_2[i]].remove(i)
            position_dict_2[arm_1[i]].append(i)
            arm_2[i] = arm_1[i]
        
        else: 
            position_dict_1[arm_1[i]].remove(i)
            position_dict_1[arm_2[i]].append(i)
            arm_1[i] = arm_2[i]

         
   
     
def simulate_one_generation(arm_1, arm_2, conversion_rate, substitution_rate_AT, substitution_rate_CG, position_dict_arm_1, position_dict_arm_2):
    
    #Select positions to mutate
    #You have n = 37000, p = mutation rate; draw once from binomial distribution to see how many mutations you get (will usually be 0)
    number_mut_A_1 = binom.rvs(len(position_dict_arm_1["A"]), substitution_rate_AT, size = 1)
    number_mut_T_1 = binom.rvs(len(position_dict_arm_1["T"]), substitution_rate_AT, size = 1)
    number_mut_C_1 = binom.rvs(len(position_dict_arm_1["C"]), substitution_rate_CG, size = 1)
    number_mut_G_1 = binom.rvs(len(position_dict_arm_1["G"]), substitution_rate_CG, size = 1)
    
    number_mut_A_2 = binom.rvs(len(position_dict_arm_2["A"]), substitution_rate_AT, size = 1)
    number_mut_T_2 = binom.rvs(len(position_dict_arm_2["T"]), substitution_rate_AT, size = 1)
    number_mut_C_2 = binom.rvs(len(position_dict_arm_2["C"]), substitution_rate_CG, size = 1)
    number_mut_G_2 = binom.rvs(len(position_dict_arm_2["G"]), substitution_rate_CG, size = 1)
    
    
    dictionary_1 = {"A": [position_dict_arm_1["A"], number_mut_A_1], "T": [position_dict_arm_1["T"], number_mut_T_1],
    "C": [position_dict_arm_1["C"], number_mut_C_1], "G": [position_dict_arm_1["G"], number_mut_G_1]}
    
    dictionary_2 = {"A": [position_dict_arm_2["A"], number_mut_A_2], "T": [position_dict_arm_2["T"], number_mut_T_2],
    "C": [position_dict_arm_2["C"], number_mut_C_2], "G": [position_dict_arm_2["G"], number_mut_G_2]}
    
    for base in dictionary_1.keys():
        
        possible_positions = dictionary_1[base][0]
        number_mut = dictionary_1[base][1]
        
        if dictionary_1[base][1]!=0:
        
            mutate_arm(arm_1, possible_positions, number_mut, position_dict_arm_1)
    
    for base in dictionary_2.keys():
        
        possible_positions = dictionary_2[base][0]
        number_mut = dictionary_2[base][1]
        
        if dictionary_2[base][1]!=0:
        
            mutate_arm(arm_2, possible_positions, number_mut, position_dict_arm_2)
    

    #Select positions for gene conversion
    number_conv = binom.rvs(len(arm_1), conversion_rate, size = 1)
    
    if number_conv != 0:
        
        convert_arm(arm_1, arm_2, number_conv, position_dict_arm_1, position_dict_arm_2)
        
        
"""Main script"""

#Get mutation spectrum
substitution_matrix = generate_substitution_matrix(relative_frequencies)
substitution_rate_branch34_AT, substitution_rate_branch34_CG = calculate_substitution_rates_AT_CG(relative_frequencies, substitution_rate_branch34)
substitution_rate_branch1_AT, substitution_rate_branch1_CG = calculate_substitution_rates_AT_CG(relative_frequencies, substitution_rate_branch1)
substitution_rate_branch2_AT, substitution_rate_branch2_CG = calculate_substitution_rates_AT_CG(relative_frequencies, substitution_rate_branch2)

#Initialize arms
bases = ["A","T","G","C"]
p = [AT_content/2, AT_content/2, GC_content/2, GC_content/2]

np.random.seed(iteration) 
initial_arm_1 = np.random.choice(bases, size=length, p=p)

#Make Arm 2 different 
initial_arm_2 = []

p_same = 1-fraction_differences
p_diff = fraction_differences/3

difference_probs = {"A": [p_same, p_diff, p_diff, p_diff], 
"T": [p_diff, p_same, p_diff, p_diff], 
"G": [p_diff, p_diff, p_same, p_diff], 
"C": [p_diff, p_diff, p_diff, p_same]}

for i in initial_arm_1:
    
    choice = np.random.choice(bases, p=difference_probs[i])
    initial_arm_2.append(choice)
      
#Make a separate identical variables for human-chimp arms, and rhesus arms (can't link by variable name, or else they will update together)
arm_1_human_chimp = [i for i in initial_arm_1]     
arm_2_human_chimp = [i for i in initial_arm_2]  
arm_1_rhesus = [i for i in initial_arm_1]
arm_2_rhesus = [i for i in initial_arm_2]
arm_1_ancestral = [i for i in initial_arm_1]
arm_2_ancestral = [i for i in initial_arm_2]  

#Make lists of positions for each base
position_dict_rhesus_arm_1 = separate_sequence_by_base(initial_arm_1)
position_dict_rhesus_arm_2 = separate_sequence_by_base(initial_arm_2)

position_dict_human_chimp_arm_1 = separate_sequence_by_base(initial_arm_1)
position_dict_human_chimp_arm_2 = separate_sequence_by_base(initial_arm_2)


#Rhesus simulations:  1450000 generations (1 generation / 20 years * 29 million years)
for i in list(range(0,1450000 )):
    
    simulate_one_generation(arm_1_rhesus, arm_2_rhesus, conversion_rate_branch1, substitution_rate_branch1_AT, 
    substitution_rate_branch1_CG, position_dict_rhesus_arm_1, position_dict_rhesus_arm_2)

#Human-chimpanzee common ancestor simulations:  1100000 generations (1 generation / 20 years * 22 million years)
for i in list(range(0,1100000)):
    
    simulate_one_generation(arm_1_human_chimp, arm_2_human_chimp, conversion_rate_branch2, substitution_rate_branch2_AT, 
    substitution_rate_branch2_CG, position_dict_human_chimp_arm_1, position_dict_human_chimp_arm_2)
    

#Update your arms:  Chimp and human are now diverged!
arm_1_human = [i for i in arm_1_human_chimp]  
arm_2_human = [i for i in arm_2_human_chimp]  
arm_1_chimp = [i for i in arm_1_human_chimp]  
arm_2_chimp = [i for i in arm_2_human_chimp]  

#Update your position dictionaries
position_dict_human_arm_1 = separate_sequence_by_base(arm_1_human)
position_dict_human_arm_2 = separate_sequence_by_base(arm_2_human)

position_dict_chimp_arm_1 = separate_sequence_by_base(arm_1_chimp)
position_dict_chimp_arm_2 = separate_sequence_by_base(arm_2_chimp)


#Human and chimpanzee simulations:  350000 (1 generation / 20 years * 7 million years)
for i in list(range(0,350000)):
    
    simulate_one_generation(arm_1_human, arm_2_human, conversion_rate_branch34, substitution_rate_branch34_AT,
    substitution_rate_branch34_CG, position_dict_human_arm_1, position_dict_human_arm_2)

for i in list(range(0,350000)):
    
    simulate_one_generation(arm_1_chimp, arm_2_chimp, conversion_rate_branch34, substitution_rate_branch34_AT,
    substitution_rate_branch34_CG, position_dict_chimp_arm_1, position_dict_chimp_arm_2)
      
#Write to files
arm_1_human = "".join(arm_1_human)
arm_2_human = "".join(arm_2_human)

arm_1_chimp = "".join(arm_1_chimp)
arm_2_chimp = "".join(arm_2_chimp)

arm_1_rhesus = "".join(arm_1_rhesus)
arm_2_rhesus = "".join(arm_2_rhesus)

arm_1_ancestral = "".join(arm_1_ancestral)
arm_2_ancestral = "".join(arm_2_ancestral)

if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

with open("{0}/Iteration_{1}.fasta".format(output_dir, iteration),"w") as new:
    new.write(""">Arm_1_rhesus\n{0}\n>Arm_2_rhesus\n{1}\n>Arm_1_human\n{2}\n>Arm_2_human\n{3}\n>Arm_1_chimp\n
    {4}\n>Arm_2_chimp\n{5}\n>Arm_1_ancestral\n{6}\n>Arm_2_ancestral\n{7}""".format(arm_1_rhesus,arm_2_rhesus, 
    arm_1_human, arm_2_human, arm_1_chimp, arm_2_chimp, arm_1_ancestral, arm_2_ancestral)) 
