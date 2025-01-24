# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 16:23:40 2020

@author: Saminur Islam
"""

import os
import csv
import numpy as np
saving_directory='relation.csv'
directory_path='Input Data\\maxdist_str'
drugbank_pdb_dir = 'drugTable.csv'

dictionary_chains= dict()
for file in os.listdir(directory_path):
    if file.endswith(".csv"):
        path=os.path.join(directory_path, file)
        with open(path) as csvfile:
            readCSV = csv.reader(csvfile)
            for row in readCSV:
                pdbchain=row[0]
                pdbid=pdbchain[:(len(pdbchain)-1)]
                chainId = str(row[0][-1]).upper()
                dna_seq=row[2:]
                dna_seq_str=''.join(dna_seq)
                if dictionary_chains.get(pdbid):
                    listVal=dictionary_chains.get(pdbid)
                    if (chainId,dna_seq_str) in listVal: continue
                    
                    if type(listVal)==tuple:
                        plusVal=[listVal]+[(chainId,dna_seq_str)]
                    else:
                        plusVal=listVal+[(chainId,dna_seq_str)]
    
                    dictionary_chains[pdbid]=plusVal
                else:
                    dictionary_chains[pdbid]=(chainId,dna_seq_str)
                    

# pdbid_list= []

# for ii in dictionary_chains:
#     pdbid_list.append(ii)
    
# new_pdb_list = []
# drugid_list = []
# dictionary = dict()
# with open(drugbank_pdb_dir,'r') as drugbankFile:
#     count = 0
#     readCSV = csv.reader(drugbankFile)
#     for row in readCSV:
#         if count == 0: 
#             count=count+1
#             continue
#         if len(row) < 6: continue 
#         pdb1=row[6]
#         if len(row) < 8: 
#             pdb2=''
             
#         else:
#             pdb2=row[8]
#         if len(row) < 10: 
#             pdb3=''
             
#         else:
#             pdb3=row[10]    
#         DrugId = row[2]  
#         drugid_list.append(DrugId)
#         if dictionary.get(DrugId):
#             listVal=dictionary.get(DrugId)
#             if type(listVal)==str:
#                 if pdb1 in pdbid_list:
#                     plusVal=[listVal]+[pdb1]
#                     new_pdb_list.append(pdb1)
#                 elif pdb2 in pdbid_list:
#                     plusVal=[listVal]+[pdb2]
#                     new_pdb_list.append(pdb2)
#                 elif pdb3 in pdbid_list:
#                     plusVal=[listVal]+[pdb3]
#                     new_pdb_list.append(pdb3)
                
#             else:
#                 if pdb1 in pdbid_list:
#                     plusVal=listVal+[pdb1]
#                     new_pdb_list.append(pdb1)
#                 elif pdb2 in pdbid_list:
#                     plusVal=listVal+[pdb2]
#                     new_pdb_list.append(pdb2)
#                 elif pdb3 in pdbid_list:
#                     plusVal=listVal+[pdb3]   
#                     new_pdb_list.append(pdb3)
#             removeDupl=list(dict.fromkeys(plusVal))
#             dictionary[DrugId]=removeDupl
#         else:
#             if pdb1 in pdbid_list:
#                 dictionary[DrugId]=pdb1
#                 new_pdb_list.append(pdb1)
#             elif pdb2 in pdbid_list:
#                 dictionary[DrugId]=pdb2
#                 new_pdb_list.append(pdb2)
#             elif pdb3 in pdbid_list:
#                 dictionary[DrugId]=pdb3
#                 new_pdb_list.append(pdb3)

import uniprot_pdb_mapping as pdbdbr
dictionary = pdbdbr.pdbid_drugid_relation()

# new_pdb_list=list(dict.fromkeys(new_pdb_list))
drugid_list=list(dictionary.keys())
# print(drugid_list.index('DB00098'))
# drugid_list = drugid_list[106:]
similarity_directory='drug_pdb_similarity_modified.csv'
columns = ['drug id One','drug id two','Levenshtein distance min','Levenshtein distancee max','Levenshtein distance avg', 'Levenshtein distance weighted avg',
           'Levenshtein distance min avg','Levenshtein distancee max avg']

import textdistance
import math

def find_the_similarites(pdbFirstSeqInfo,pdbSecondSeqInfo):
    max_simi = 0
    min_simi = 1
    avg_simi = 0
    exponetial_avg = 0
    list_val = []
    exp_list = []
    for ii in range(len(pdbFirstSeqInfo)):
        if type(pdbFirstSeqInfo)==tuple:
            (X1,Y1) = pdbFirstSeqInfo
        else:
            (X1,Y1) = pdbFirstSeqInfo[ii]
        for jj in range(len(pdbSecondSeqInfo)):
            if type(pdbSecondSeqInfo)==tuple:
                (X2,Y2) = pdbSecondSeqInfo
            else:
                (X2,Y2) = pdbSecondSeqInfo[jj] 
            levn = textdistance.levenshtein.normalized_similarity(Y1, Y2)
            list_val.append(levn)
            if levn >  max_simi:
                max_simi = levn
                
            if levn < min_simi:
                min_simi = levn
            avg_simi = avg_simi + levn
            try:
                ans = math.exp(levn)
            except OverflowError:
                ans = float('inf')
            exp_list.append(ans)
            try:
                ans2 = math.exp(ans)
            except OverflowError:
                ans2 = float('inf')
            exponetial_avg = exponetial_avg + ans2
    avg_simi = avg_simi/ (len(pdbFirstSeqInfo) * len(pdbSecondSeqInfo))
    weight_list  = np.divide(exp_list, exponetial_avg)
    final_weighted_avg = 0
    # print(len(weight_list))
    # print(len(list_val))
    for i in range(len(list_val)):
        final_weighted_avg = final_weighted_avg + list_val[i] * weight_list[i]
    
    calulted_value_list = [min_simi,max_simi,avg_simi,final_weighted_avg/len(list_val)]
    return calulted_value_list

def calculate_final_similarity_values(list_val):
    max_simi = 0
    min_simi = 1
    avg_simi = 0
    exponetial_avg = 0
    max_avg = 0
    min_avg = 0
    exp_list = []
    for i in range(len(list_val)):
        mins= list_val[i][0] 
        maxs= list_val[i][1] 
        avgs = list_val[i][2] 
        weighted_avgs = list_val[i][3]
        max_avg = max_avg +  maxs
        min_avg = min_avg + mins
        if max_simi < maxs:
            max_simi= maxs
        if min_simi > mins:
            min_simi = mins
        avg_simi = avg_simi + avgs           
        exponetial_avg = exponetial_avg + weighted_avgs
    
    if len(list_val) == 0: return [0,0,0,0,0,0]
    
    avg_simi = avg_simi / len(list_val)
    max_avg = max_avg / len(list_val)
    min_avg = min_avg / len(list_val)
    exponetial_avg = exponetial_avg/ len(list_val)
    return [min_simi,max_simi,avg_simi,exponetial_avg, min_avg,max_avg]
        
with open(similarity_directory, 'w',newline='') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        
        for i in range(len(drugid_list)):
            firstBool = False
            if  drugid_list[i] in dictionary:
                pdb_list = dictionary.get(drugid_list[i])
            else:
                continue
            
            if type(pdb_list)==str:
                firstBool = True
            for j in range(i+1,len(drugid_list)):  
                secndBool = False
                if  drugid_list[j] in dictionary:
                     pdb_list_scnd = dictionary.get(drugid_list[j])
                else:
                    continue
                
                if type(pdb_list_scnd)==str:
                    secndBool = True
                    
                if firstBool == secndBool == True: 
                    pdbFirstSeqInfo = dictionary_chains.get(pdb_list)
                    pdbSecondSeqInfo = dictionary_chains.get(pdb_list_scnd)
                    if pdbFirstSeqInfo is None or pdbSecondSeqInfo is None: continue
                    calculate_similarities = find_the_similarites(pdbFirstSeqInfo,pdbSecondSeqInfo)
                    row = [drugid_list[i],drugid_list[j]] + calculate_similarities+[calculate_similarities[0]]+[calculate_similarities[1]]
                
                elif firstBool== True and secndBool == False:
                    pdbFirstSeqInfo = dictionary_chains.get(pdb_list)
                    if pdbFirstSeqInfo is None: continue
                    list_val = []
                    for mm in range(len(pdb_list_scnd)):
                        if dictionary_chains.get(pdb_list_scnd[mm]) is None: continue
                        if pdb_list_scnd[mm] in dictionary_chains:
                            pdbSecondSeqInfo = dictionary_chains.get(pdb_list_scnd[mm])
                        else: 
                            continue
                        calculate_similarities = find_the_similarites(pdbFirstSeqInfo,pdbSecondSeqInfo)
                        list_val.append(calculate_similarities)
                        if len(list_val) ==0: 
                            continue
                    processed = calculate_final_similarity_values(list_val)
                    row = [drugid_list[i],drugid_list[j]] + processed
                                                                  
                elif firstBool== False and secndBool == True:
                    pdbSecondSeqInfo = dictionary_chains.get(pdb_list_scnd)
                    if pdbSecondSeqInfo is None: continue
                    list_val = []
                    for mm in range(len(pdb_list)):
                        if dictionary_chains.get(pdb_list[mm]) is None: continue
                        pdbFirstSeqInfo = dictionary_chains.get(pdb_list[mm])
                        calculate_similarities = find_the_similarites(pdbFirstSeqInfo,pdbSecondSeqInfo)
                        list_val.append(calculate_similarities)
                        if len(list_val) ==0: 
                            continue
                    processed = calculate_final_similarity_values(list_val)
                    row = [drugid_list[i],drugid_list[j]] + processed
                else:
                    list_val = []
                    for mm in range(len(pdb_list)):
                        if dictionary_chains.get(pdb_list[mm]) is None: continue
                        if pdb_list[mm] in dictionary_chains:
                            pdbFirstSeqInfo = dictionary_chains.get(pdb_list[mm])
                        for nn in range(len(pdb_list_scnd)):
                            if dictionary_chains.get(pdb_list_scnd[nn]) is None: continue
                            if pdb_list_scnd[nn] in dictionary_chains:
                                pdbSecondSeqInfo = dictionary_chains.get(pdb_list_scnd[nn])
                            else: 
                                continue
                            calculate_similarities = find_the_similarites(pdbFirstSeqInfo,pdbSecondSeqInfo)
                            list_val.append(calculate_similarities)
                            if len(list_val) ==0: 
                                continue
                    processed = calculate_final_similarity_values(list_val)
                    row = [drugid_list[i],drugid_list[j]] + processed
                    print(row)
                    writer.writerow(row)
# import pandas as pd
# import textdistance
# similarity_directory='G:\\Study\\Research WVU\\PDB\\similarity_drugbank\\'
# for i in range(len(new_pdb_list)):
# #    print(drugList[i])
#     columns = ['pdbId one','PdbId two','similarity','chain one','chain two']
#     file = new_pdb_list[i]+'.csv'
#     path=os.path.join(similarity_directory, file)
#     with open(path, 'w',newline='') as f:
#         writer = csv.writer(f)
#         writer.writerow(columns)
#         for j in range(i+1,len(new_pdb_list)):
# #        print(protientargetList[j])

#             pdbFirstSeqInfo = dictionary_chains.get(new_pdb_list[i])
#             pdbSecondSeqInfo = dictionary_chains.get(new_pdb_list[j])
#             # print(new_pdb_list[i])
#             # print(pdbFirstSeqInfo)
#             max_simi=0
#             C1=''
#             C2=''
#             for ii in range(len(pdbFirstSeqInfo)):
#                 if type(pdbFirstSeqInfo)==tuple:
#                     (X1,Y1) = pdbFirstSeqInfo
#                 else:
#                     (X1,Y1) = pdbFirstSeqInfo[ii]
#                 for jj in range(len(pdbSecondSeqInfo)):
#                     if type(pdbSecondSeqInfo)==tuple:
#                         (X2,Y2) = pdbSecondSeqInfo
#                     else:
#                         (X2,Y2) = pdbSecondSeqInfo[jj] 
#                     # if len(pdbSecondSeqInfo[jj])!=2:
#                     #     # print(pdbid_list[j])
#                     #     continue
#                     # (X2,Y2) = pdbSecondSeqInfo[jj]
#                     levn = textdistance.levenshtein.normalized_similarity(Y1, Y2)
#                     if levn >  max_simi:
#                         max_simi = levn
#                         C1 = X1
#                         C2 = X2
#             row = [new_pdb_list[i],new_pdb_list[j],max_simi,C1,C2]
#             writer.writerow(row)
#     print("complete file name: "+str(file))
            
# with open(saving_directory, 'w',newline='') as f:
#     writer = csv.writer(f)
#     for key in dictionary.keys():
#         row=[key]+dictionary[key]
#         writer.writerow(row)
        

  
        
