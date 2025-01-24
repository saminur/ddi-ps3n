# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 18:06:15 2020

@author: Saminur Islam
"""


import os
import csv
import return_pdb_from_drugbank as pdbdb
import math
import numpy as np
import uniprot_pdb_mapping as pdbdbr


saving_directory='relation.csv'
directory_path='kmers_directory'
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
                kmers_list = []
                if len(row) < 2: continue
            
                for i in range(2,len(row),2):
                    k_tuple = (row[i],int(row[i+1]))
                    kmers_list.append(k_tuple)
                    
                dna_seq=row[2:]
                dna_seq_str=''.join(dna_seq)
                if dictionary_chains.get(pdbid):
                    listVal=dictionary_chains.get(pdbid)
                    # if (chainId,dna_seq_str) in listVal: continue
                    plusVal=listVal+kmers_list
                    d = {x:0 for x, _ in plusVal} 
                    for name, num in plusVal: d[name] += num 
                    Output = list(map(tuple, d.items())) 
                    dictionary_chains[pdbid]=Output
                else:
                    dictionary_chains[pdbid]=kmers_list




dictionary = pdbdbr.function_uniprot_drug_pdb_relation()
new_list = list(dictionary.keys())

import similarity_kmer_profiles as skp

skp.calculating_profile_similarity(dictionary,dictionary_chains,list(dictionary.keys()))
# print(dictionary_chains['4XPH'])
# print(dictionary_chains['5CXV'])
# a = [1, 2, 3]
# b = [1, 2, 3, 4]
# from itertools import zip_longest
# c= [x + y for x, y in zip_longest(a, b, fillvalue=0)]
# print(c)
# print(c*.5)
# firstkmer=dictionary_chains['4LLR']
# secondkmer= dictionary_chains['4XCS']

# sortedOne = sorted(firstkmer, key=lambda firstkmer: firstkmer[0])
# sortedTwo = sorted(secondkmer, key=lambda secondkmer: secondkmer[0])
# mergedList = sortedOne + sortedTwo
# d = {x:0 for x, _ in mergedList} 
# d_first ={x:0 for x, _ in sortedOne} 
# for name, num in sortedOne: d_first[name] += num 
# d_scnd ={x:0 for x, _ in sortedTwo} 
# for name, num in sortedTwo: d_scnd[name] += num 
# for name, num in mergedList: d[name] += num 
# Output = list(map(tuple, d.items()))

# convertedList1 =[]
# convertedList2 =[]
# for ii in range(len(Output)):
#     if Output[ii][0] in d_first:
#         val = d_first.get(Output[ii][0])
#         frequency = float(val) / float(Output[ii][1])
#         convertedList1.append(frequency)
#     else:
#         convertedList1.append(0.0)
    
#     if Output[ii][0] in d_scnd:
#         val2 = d_scnd.get(Output[ii][0])
#         frequency2 = float(val2) / float(Output[ii][1])
#         convertedList2.append(frequency2)
#     else:
#         convertedList2.append(0.0)
        

# print(1.0/euclidean_distance(convertedList1,convertedList2))
# print(cosine_similarity(convertedList1,convertedList2))
# values1 = np.asarray(convertedList1)
# values2 = np.asarray(convertedList2)
# # print(1/KL(values1,values2))


# print(dictionary_chains['4LLR'])
# print(dictionary_chains['4XCS'])
# list1 = [('a', 1), ('b', 0), ('c', 0)]
# list2 = [('a', 5), ('c', 3)]
# d1 = dict(list1)
# d1.update(dict(list2))
# print(list(d1.items()))
