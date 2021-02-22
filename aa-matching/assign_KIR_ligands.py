# This is a script to assign KIR ligand categories 

import pandas as pd
import re
import csv
import os
import sys
from collections import defaultdict

from aa_matching_msf import *
aa_mm = AAMatch(dbversion=3420)

input_file = "impute_cases.AFA.csv"

# set filepath for loading input files
os.chdir(os.path.dirname(os.path.abspath(__file__)))
path = os.getcwd()

# load the input file specified in the command line and choose the relevant columns
impute_cases_POP = pd.read_csv(os.path.join(path, input_file), header=None)

# assign column names to impute_cases_AFA
impute_cases_POP.columns = ['ID', 'Imputation','Haplotype_1','Haplotype_2','Frequency']

# expand haplotype data to include columns for each allele
impute_cases_POP[['A1','C1','B1','DRBX_1','DRB1_1','DQA1_1','DQB1_1','DPA1_1','DPB1_1']] = impute_cases_POP.Haplotype_1.str.split("~",expand=True)
impute_cases_POP[['A2','C2','B2','DRBX_2','DRB1_2','DQA1_2','DQB1_2','DPA1_2','DPB1_2']] = impute_cases_POP.Haplotype_2.str.split("~",expand=True)

# call the dataframe defined above so you can append the KIR ligand assignments
final_df = pd.DataFrame(impute_cases_POP, columns = ['ID', 'Imputation','Haplotype_1','Haplotype_2','Frequency','A1','C1','B1','DRBX_1','DRB1_1','DQA1_1','DQB1_1','DPA1_1','DPB1_1','A2','C2','B2','DRBX_2','DRB1_2','DQA1_2','DQB1_2','DPA1_2','DPB1_2'])

# KIR ligand assignments for HLA-A allele on column "A1"

result_KIR_A_Ligand1 = []
for i in range(len(impute_cases_POP)):
    allele1 = (impute_cases_POP.loc[i, "A1"])
    position1 = 83
    allele_substring_1 = "A*11:"
    allele_substring_2 = "A*03:"
    AA_83 = aa_mm.getAAposition(allele1,position1)
    if allele_substring_1 in allele1:
        result_KIR_A_Ligand1.append("KIR_Ligand_A" + "_A311")
    elif allele_substring_2 in allele1:
        result_KIR_A_Ligand1.append("KIR_Ligand_A" + "_A311")
    elif AA_83 == "R":
        result_KIR_A_Ligand1.append("KIR_Ligand_A" + "_Bw4")
    else:
        result_KIR_A_Ligand1.append("KIR_Ligand_A" + "_X")

# KIR ligand assignments for HLA-A allele on column "A2"

result_KIR_A_Ligand2 = []
for i in range(len(impute_cases_POP)):
    allele1 = (impute_cases_POP.loc[i, "A2"])
    position1 = 83
    allele_substring_1 = "A*11:"
    allele_substring_2 = "A*03:"
    AA_83 = aa_mm.getAAposition(allele1,position1)
    if allele_substring_1 in allele1:
        result_KIR_A_Ligand2.append("KIR_Ligand_A" + "_A311")
    elif allele_substring_2 in allele1:
        result_KIR_A_Ligand2.append("KIR_Ligand_A" + "_A311")
    elif AA_83 == "R":
        result_KIR_A_Ligand2.append("KIR_Ligand_A" + "_Bw4")
    else:
        result_KIR_A_Ligand2.append("KIR_Ligand_A" + "_X")

# KIR ligand assignments for HLA-C allele on column "C1"

result_KIR_C_Ligand1 = []
for i in range(len(impute_cases_POP)):
    allele1 = (impute_cases_POP.loc[i, "C1"])
    position1 = 80
    position2 = 76
    AA_80 = aa_mm.getAAposition(allele1,position1)
    AA_76 = aa_mm.getAAposition(allele1,position2)
    print(allele1)
    print(AA_80)
    
    if AA_80 == "N" and AA_76 == "V":
        result_KIR_C_Ligand1.append("KIR_Ligand_C" + "_C1")
    elif AA_80 == "K" and AA_76 == "V":
        result_KIR_C_Ligand1.append("KIR_Ligand_C" + "_C2")
    else:
        result_KIR_C_Ligand1.append("KIR_Ligand_C" + "_X")

# KIR ligand assignments for HLA-C allele on column "C2"

result_KIR_C_Ligand2 = []
for i in range(len(impute_cases_POP)):
    allele1 = (impute_cases_POP.loc[i, "C2"])
    position1 = 80
    position2 = 76
    AA_80 = aa_mm.getAAposition(allele1,position1)
    AA_76 = aa_mm.getAAposition(allele1,position2)
    if AA_80 == "N" and AA_76 == "V":
        result_KIR_C_Ligand2.append("KIR_Ligand_C" + "_C1")
    elif AA_80 == "K" and AA_76 == "V":
        result_KIR_C_Ligand2.append("KIR_Ligand_C" + "_C2")
    else:
        result_KIR_C_Ligand2.append("KIR_Ligand_C" + "_X")

# KIR ligand assignments for HLA-B allele on column "B1"

result_KIR_B_Ligand1 = []
for i in range(len(impute_cases_POP)):
    allele1 = (impute_cases_POP.loc[i, "B1"])
    position1 = 80
    position2 = 76
    position3 = 83
    AA_80 = aa_mm.getAAposition(allele1,position1)
    AA_76 = aa_mm.getAAposition(allele1,position2)
    AA_83 = aa_mm.getAAposition(allele1,position3)
    if AA_80 == "N" and AA_76 == "V":
        result_KIR_B_Ligand1.append("KIR_Ligand_B" + "_C1")
    elif AA_83 == "R":
        result_KIR_B_Ligand1.append("KIR_Ligand_B" + "_Bw4")
    else:
        result_KIR_B_Ligand1.append("KIR_Ligand_B" + "_X")

# KIR ligand assignments for HLA-B allele on column "B2"

result_KIR_B_Ligand2 = []
for i in range(len(impute_cases_POP)):
    allele1 = (impute_cases_POP.loc[i, "B2"])
    position1 = 80
    position2 = 76
    position3 = 83
    AA_80 = aa_mm.getAAposition(allele1,position1)
    AA_76 = aa_mm.getAAposition(allele1,position2)
    AA_83 = aa_mm.getAAposition(allele1,position3)
    if AA_80 == "N" and AA_76 == "V":
        result_KIR_B_Ligand2.append("KIR_Ligand_B" + "_C1")
    elif AA_83 == "R":
        result_KIR_B_Ligand2.append("KIR_Ligand_B" + "_Bw4")
    else:
        result_KIR_B_Ligand2.append("KIR_Ligand_B" + "_X")

# add the KIR ligand categories generated by the loops to the original dataframe

final_df["KIR_A_Ligand1"] = result_KIR_A_Ligand1
final_df["KIR_A_Ligand2"] = result_KIR_A_Ligand2
final_df["KIR_C_Ligand1"] = result_KIR_C_Ligand1
final_df["KIR_C_Ligand2"] = result_KIR_C_Ligand2
final_df["KIR_B_Ligand1"] = result_KIR_B_Ligand1
final_df["KIR_B_Ligand2"] = result_KIR_B_Ligand2

print(final_df)

# output final_df to csv based on population name in current working directory
# cwd = os.getcwd()
# output_path = cwd + "/KIR_ligand_" + input_file
# final_df.to_csv(output_path, sep = ",", index = False)
