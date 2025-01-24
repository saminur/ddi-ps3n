# -*- coding: utf-8 -*-
"""
Created on Wed May 20 15:09:37 2020

@author: Saminur Islam
"""

from Bio import pairwise2
from Bio.pairwise2 import format_alignment 
from collections import Counter
from math import sqrt
# from similarity.jarowinkler import JaroWinkler
from Bio import SeqIO
import csv
import kmers_jaccard as kmers
from pandas import read_excel
import pandas as pd
#import numpy as np
#from Bio.Blast import NCBIWWW
fasta_file="protein.fasta"

from pandas import read_excel

def stringTovector(sequence):
    
    # count the characters in word
    counts = Counter(sequence)
    # precomputes a set of the different characters
    soc = set(counts)
    # precomputes the "length" of the word vector
    lengthSeq = sqrt(sum(c*c for c in counts.values()))
    # return a tuple
    return counts, soc, lengthSeq

def cosine_similarity(str1, str2):
    v1=stringTovector(str1)
    v2=stringTovector(str2)
    # which characters are common to the two words
    common = v1[1].intersection(v2[1])
    # by definition of cosine distance we have
    return sum(v1[0][ch]*v2[0][ch] for ch in common)/v1[2]/v2[2]

my_sheet = 'Sheet1' # change it to your sheet name
file_name = 'pathways.xlsx'  # change it to the name of your excel file
df = read_excel(file_name, sheet_name = my_sheet)
dictionary=dict()
drugList=[]
protientargetList=[]
for i in range(len(df)):
    DrugIds=str(df.iloc[i][3]).split(",")
    uniprotIds=str(df.iloc[i][4]).split(",")
    drugList=drugList+DrugIds
    protientargetList=protientargetList+uniprotIds
    for k in range(len(DrugIds)):
        did=DrugIds[k]
        if dictionary.get(did):
            listVal=dictionary.get(did)
            plusVal=listVal+uniprotIds
            removeDupl=list(dict.fromkeys(plusVal))
            dictionary[did]=removeDupl
        else:
            dictionary[did]=uniprotIds

drugList=list(dict.fromkeys(drugList))
print(drugList.index('DB02083'))
# drugList = drugList[40:]
protientargetList=list(dict.fromkeys(protientargetList))          
seq_record = next(SeqIO.parse(open(fasta_file),'fasta'))
print(seq_record)
print(seq_record.seq)
seq_obj_list={}
list_of_index=[]
with open(fasta_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        stringformat=str(record.id)
        prot_id=stringformat.split("|")
        list_of_index.append(prot_id[1])
        seq_obj_list[prot_id[1]]=str(record.seq)
        

#k-mer profile similarity for cosine, euclidean and JS divergence

kmer_protein_seq_dic =dict()
for key in seq_obj_list:
    suffix_arr = kmers.computeSuffixArray(seq_obj_list.get(key))
    kmer_touple_list = kmers.findValidPairs(suffix_arr,len(seq_obj_list.get(key)),3)  
    kmer_protein_seq_dic[key] = kmer_touple_list


import similarity_functions_proteins as skp    
skp.calculating_profile_similarity(dictionary,kmer_protein_seq_dic,drugList)

# gg = seq_obj_list.get('P19113')
# print(gg)
# suffix_arr = kmers.computeSuffixArray(gg)
# kmer_touple_list = kmers.findValidPairs(suffix_arr,len(gg),6)    
# print(kmer_touple_list)


 

# another_path="H:\\Study\\Spring2020Semester\\Research\\New DrugnMetabolism\\DrugBank\\drugbank_approved_target_polypeptide_sequences.fasta\\protein_similarity_cosine.csv"
# columns_an=['protein0','protien1','similarity_score']
# with open(another_path, 'w',newline='') as f:
#     writer = csv.writer(f)
#     writer.writerow(columns_an)
#     for ii in range(len(list_of_index)):
#         first_seq=seq_obj_list.get(list_of_index[ii])
#         for jj in range(ii+1,len(list_of_index)):
#             scnd_seq=seq_obj_list.get(list_of_index[jj])
#             similarity_value = cosine_similarity(first_seq, scnd_seq)
#             row=[list_of_index[ii],list_of_index[jj],similarity_value]
#             writer.writerow(row)
#            df_final.at[list_of_index[ii],list_of_index[jj]]=similarity_value
            


#another_path="H:\\Study\\Spring2020Semester\\Research\\New DrugnMetabolism\\DrugBank\\drugbank_approved_target_polypeptide_sequences.fasta\\protein_similarity_ano.csv"
#columns_an=['drug id0', 'drug id1', 'similarity_score']
#with open(another_path, 'w',newline='') as f:
#    writer = csv.writer(f)
#    writer.writerow(columns_an)
#    for (id0, fp0), (id1, fp1) in itertools.combinations(dictionary.items(), 2):
#        sm_scor=similarity_score(fp0,fp1)
#        if sm_scor==0.0: continue
#        row=[id0,id1,sm_scor]
#        writer.writerow(row)
