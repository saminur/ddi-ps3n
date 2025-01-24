# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 14:22:08 2020

@author: Saminur Islam
"""

import os
import csv
import numpy as np
import pandas as pd

p1 ='ddi_kmers.csv'
df = pd.read_csv(p1)
df = df.set_index('Unnamed: 0')
d1 = list(df.keys())
for i in range(len(d1)):
    df.at[d1[i],d1[i]] = 1

df.to_csv('ddi_kmers.csv')
t1 = 'data\\'

saving_dir = 'SNF matrix\\'

feature_path = 'drug_drug_int.csv'
drug_drug_int = pd.read_csv(feature_path)
drug_drug_int = drug_drug_int.set_index('Unnamed: 0')
drugList_int = list(drug_drug_int.keys())
dir_path = 'dataset\\'
ddmain = set(drugList_int)
pdb_kmer_path = dir_path + 'pdb_similarity_deepddi_new.csv'
df_pdb_kmer = pd.read_csv(pdb_kmer_path)
drugList_pdb = list(dict.fromkeys(df_pdb_kmer['drug id One']))
drugList = list(ddmain.intersection(drugList_pdb))

drug_drug_int['DB00001']['DB00002']

df_drug  = pd.DataFrame()
for i in drugList:
    for j in drugList:
        #print(drug_drug_int[i][j])
        df_drug.at[i,j] = drug_drug_int[i][j]
df_drug.to_csv("ddi_structure.csv")

for filename in os.listdir(t1):
    if "ddi" in filename:
        continue
    else:
        path = os.path.join(t1,filename)
        s_dir = os.path.join(saving_dir,filename)
        mat_ = []
        df = pd.read_csv(path)
        df = df.set_index('Unnamed: 0')
        drugList = list(df.keys())
        with open(s_dir, 'w',newline='') as f:
            writer = csv.writer(f)
            for ii in range(len(drugList)):
                for jj in range(ii+1,len(drugList)):
                    firstList = df[drugList[ii]]
                    scndList = df[drugList[jj]]
                    output = list(firstList - scndList)
                    row = [drugList[ii],drugList[jj]] + output
                    writer.writerow(row)


path = 'pairwise_similarity_indications.csv'
df = pd.read_csv(path)
df = df.set_index('Unnamed: 0')
drugList = list(df.keys())



matrices = []


t1 = 'dataset\\'
drugList = []
for filename in os.listdir(t1):
    if "ddi" in filename:
        path = os.path.join(t1,filename)
        df = pd.read_csv(path)
        df = df.set_index('Unnamed: 0')
        drugList = list(df.keys())
    else:
        path = os.path.join(t1,filename)
        mat_ = []
        with open(path) as csvfile:
            readCSV = csv.reader(csvfile)
            count = 0
            for row in readCSV:
                if count == 0:
                    count += 1
                    continue
                
                row = row[1:]
                row = ["0" if x == '' else x for x in row]
                row = list(map(np.float32,row))
                mat_.append(row)
        matrices.append(np.array(mat_))


# print(np.array(mat_).shape)  
# drugList = list(np.array(mat_).keys())                   
import snf
print(4)
affinity_networks = snf.make_affinity(matrices,metric='euclidean',K=5, mu = 0.30)
fused_network = snf.snf(affinity_networks, K=10)
print(fused_network.shape)
drugList = sorted(drugList)
import pandas as pd
df_dat = pd.DataFrame(fused_network,index=drugList,columns=drugList)
df_dat.to_csv("merged_feature.csv")



dir_ = 'pdb_cosine_avg.csv'
df = pd.read_csv(dir_)
df = df.set_index('Unnamed: 0')
import numpy as np
# define a matrix
A = df.to_numpy()
# # A = np.array([[1, 2,3], [3, 4,5], [5, 6,7]])
# print(A)
# # calculate the mean of each column
# M = np.mean(A.T, axis=1)
# print(M)
# # center columns by subtracting column means
# C = A - M
# print(C)
# # calculate covariance matrix of centered matrix
# V = np.cov(C.T)
# print(V)
# # eigendecomposition of covariance matrix
# values, vectors = np.linalg.eig(V)
# print(vectors)
# print(values)
# # project data
# P = vectors.T.dot(C.T)
# print(P.T)

from scipy.spatial.distance import pdist,squareform
k = [
     [2, 4, 7], 
     [3, 4, 7], 
     [5, 1, 3]
     ]
dm = pdist(np.array(A), 'euclidean')
dm = squareform(dm)
sum_row = np.sum(dm,axis = 1)
sum_row = np.reshape(sum_row,(-1,1))
dm = dm/sum_row

# sum_ = 0
# for i in range(len(dm)):
#     sum_ = sum_ + dm[1][i]
# sum_row_next = np.sum(dm,axis = 1)


p3 = 'DS3\\ddi_pdb.csv'
p1 ='DS1\\ddi_kmers.csv'
df = pd.read_csv(p1)
df = df.set_index('Unnamed: 0')
d1 = list(df.keys())
df = pd.read_csv(p3)
df = df.set_index('Unnamed: 0')
d3 = list(df.keys())

d4 = list(set(d1)|set(d3))
drudlist = sorted(d4)

directory = 'pharmacologically_active_target.csv'

def pdbid_drugid_relation():
    drugbank_list = []
    drugbank_uniprot = dict()
    uniprot_pdb_map = dict()
    with open(directory) as csvfile:
        readCSV = csv.reader(csvfile)
        count = 0
        for row in readCSV:
            if count == 0:
                count = count + 1
                continue;
            if row[11] == 'Humans':
                uniprot_id = row[5]
                drug_ids = str(row [12]).split(";")
            
                if type(drug_ids) == str:
                    drug_ids=drug_ids.replace(" ", "")
                    if drugbank_uniprot.get(drug_ids):
                        uniprot_list = drugbank_uniprot.get(drug_ids)
                        new_list = uniprot_list + [uniprot_id]                   
                        drugbank_uniprot[drug_ids] = new_list
                    else:
                        drugbank_uniprot[drug_ids] = [uniprot_id]
                else: 
                    for i in range(len(drug_ids)):
                        if drugbank_uniprot.get(drug_ids[i]):
                            uniprot_list = drugbank_uniprot.get(drug_ids[i])
                            if type(uniprot_list) == str:
                                new_list = [uniprot_list] + [uniprot_id]
                            else:
                                new_list = uniprot_list + [uniprot_id]
                    
                            drugbank_uniprot[drug_ids[i].replace(" ", "")] = new_list
                        else:
                            drugbank_uniprot[drug_ids[i].replace(" ", "")] = [uniprot_id]
    return drugbank_uniprot


drugbank_uniprot = pdbid_drugid_relation()     

drugIdList2 = list(drugbank_uniprot.keys())

final = list(set(drudlist) & set(drugIdList2))
from Bio import SeqIO
import kmers_jaccard as kmers
fasta_file="protein.fasta"
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
skp.calculating_profile_similarity(drugbank_uniprot,kmer_protein_seq_dic,final)

final = d3
path  = 'MergedproteinPDBSimilarity.csv'
df_proteins =  pd.read_csv(path)
df_proteins = df_proteins.set_index(['drug id One','drug id two']).T.to_dict('list')
df_rslt = pd.DataFrame()
target_ma = 'mergeSequence.csv'
final = sorted(final)
for ii in range(len(final)):
    for jj in range(len(final)):
        tuple_ = (final[ii],final[jj])
        tuple_r = (final[jj],final[ii])
        
        if final[ii] == final[jj]:
            df_rslt.at[final[ii],final[jj]] = 1   
            
        if tuple_ in df_proteins:
            df_rslt.at[final[ii],final[jj]] = df_proteins.get(tuple_)[11]
            df_rslt.at[final[jj],final[ii]] = df_proteins.get(tuple_)[11]
            
        elif tuple_r in df_proteins:
            df_rslt.at[final[jj],final[ii]] = df_proteins.get(tuple_r)[11]
            df_rslt.at[final[ii],final[jj]] = df_proteins.get(tuple_r)[11]


df_rslt.to_csv(target_ma)

feature_path = 'drug_drug_intDD.csv'
drug_drug_int = pd.read_csv(feature_path)
drug_drug_int = drug_drug_int.set_index('Unnamed: 0')

df_newDDI = pd.DataFrame()

for ii in range(len(final)):
    for jj in range(len(final)):
        if final[ii] in drug_drug_int and final[jj] in drug_drug_int:
            
            df_newDDI.at[final[ii],final[jj]] = drug_drug_int[final[ii]][final[jj]]
            df_newDDI.at[final[jj],final[ii]] = drug_drug_int[final[ii]][final[jj]]
        else: 
            df_newDDI.at[final[ii],final[jj]] = 0
            df_newDDI.at[final[jj],final[ii]] = 0




    