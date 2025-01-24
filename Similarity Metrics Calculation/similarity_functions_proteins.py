# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 23:20:02 2020

@author: Saminur Islam
"""




import csv
import math
import numpy as np
import scipy.stats
from itertools import zip_longest

similarity_directory = 'xxx_similarity_xx.csv'
columns = ['drug id One','drug id two','cosine min','cosine max','cosine avg','Cosine Weighted avg','KL divergence min','KL divergence max','KL divergence avg','KL divergence weighted avg' 
           ,'Eucledian distance min','Eucledian distance max','Eucledian distance avg','Eucledian distance weighted avg']


def euclidean_distance(x,y):
    return math.sqrt(sum(pow(a-b,2) for a, b in zip(x, y)))


def square_rooted(x):
   return round(math.sqrt(sum([a*a for a in x])),3)
  
def cosine_similarity(x,y):
 numerator = sum(a*b for a,b in zip(x,y))
 denominator = square_rooted(x)*square_rooted(y)
 return round(numerator/float(denominator),3)


""" Epsilon is used here to avoid conditional code for
checking that neither P nor Q is equal to 0. """

def KL(P,Q):
    epsilon = 0.00001
    # You may want to instead make copies to avoid changing the np arrays.
    P = P+epsilon
    Q = Q+epsilon

    divergence = np.sum(P*np.log(P/Q))
    return divergence






def jensen_shannon_distance(p, q):
    """
    method to compute the Jenson-Shannon Distance 
    between two probability distributions
    """
    p = np.array(p)
    q = np.array(q)
    # convert the vectors into numpy arrays in case that they aren't
    if len(p) > len(q):
        difference = len(p) - len(q)
        q=np.concatenate([q, np.zeros(difference)])
    else: 
        difference = len(q) - len(p)
        p=np.concatenate([p, np.zeros(difference)])


    # calculate m
    # c= [x + y for x, y in zip_longest(p, q, fillvalue=0)]
    # m = np.array(c) / 2
    m = (p + q) / 2

    # compute Jensen Shannon Divergence
    divergence = (scipy.stats.entropy(p, m) + scipy.stats.entropy(q, m)) / 2

    # compute the Jensen Shannon Distance
    distance = np.sqrt(divergence)

    return distance

def find_the_similarites(pdbFirstSeqInfo,pdbSecondSeqInfo):
    totFreq1 = 0
    totFreq2 = 0
   
    for i in range(len(pdbFirstSeqInfo)):
        totFreq1 = totFreq1+ pdbFirstSeqInfo[i][1]
    
    for j in range(len(pdbSecondSeqInfo)): 
        totFreq2 = totFreq2 + pdbSecondSeqInfo[j][1]
    
    dist_list1 = []
    dist_list2 = []
    for i in range(len(pdbFirstSeqInfo)):
        dist_list1.append(float(pdbFirstSeqInfo[i][1])/totFreq1)
        
    for i in range(len(pdbSecondSeqInfo)):
        dist_list2.append(float(pdbSecondSeqInfo[i][1])/totFreq2)
    
    euclidean_dist_sim = 1 / (euclidean_distance(dist_list1,dist_list2) + 0.00001)
    cosine_dist_sim = cosine_similarity(dist_list1,dist_list2)
    kl_divergence_sim = 1 / jensen_shannon_distance(dist_list1,dist_list2)
    return [cosine_dist_sim,kl_divergence_sim,euclidean_dist_sim]


def calculate_final_similarity_values(list_val):
    if len(list_val) == 0: 
        return [0,0,0,0]
    max_simi = list_val[0]
    min_simi = list_val[0]
    avg_simi = 0
    exponetial_avg_sum  = 0
    exp_list = []
    for i in range(len(list_val)):
        if list_val[i] > max_simi: 
            max_simi = list_val[i]
        
        if list_val[i] < min_simi:
            min_simi = list_val[i]
        avg_simi =avg_simi + list_val[i]
        try:
            ans = math.exp(list_val[i])
        except OverflowError:
            ans = float('inf')
        exp_list.append(ans)
        try:
            ans2 = math.exp(ans)
        except OverflowError:
            ans2 = float('inf')
        exponetial_avg_sum = exponetial_avg_sum + ans2
    weight_list  = np.divide(exp_list, exponetial_avg_sum)
    final_weighted_avg = 0
    # print(len(weight_list))
    # print(len(list_val))
    avg_simi = avg_simi / len(list_val)
    for i in range(len(list_val)):
        final_weighted_avg = final_weighted_avg + list_val[i] * weight_list[i]
    return [min_simi,max_simi,avg_simi,final_weighted_avg/len(list_val)]
        
        
def divide_similarity_wise(list_val):
    cosine_sim_list = []
    kl_sim_list = []
    eucl_sim_list = []
    for i in range(len(list_val)):
        cosine_sim_list.append(list_val[i][0])
        kl_sim_list.append(list_val[i][1])
        eucl_sim_list.append(list_val[i][2])
    
    cosine_final = calculate_final_similarity_values(cosine_sim_list)
    kl_final = calculate_final_similarity_values(kl_sim_list)
    # print(eucl_sim_list)
    eucl_final = calculate_final_similarity_values(eucl_sim_list)
    
    return [cosine_final,kl_final,eucl_final]

def calculating_profile_similarity(dictionary,dictionary_chains,drugid_list):
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
                    
                     row = [drugid_list[i],drugid_list[j],calculate_similarities[0],calculate_similarities[0],calculate_similarities[0],calculate_similarities[0],
                           calculate_similarities[1],calculate_similarities[1],calculate_similarities[1],calculate_similarities[1],
                           calculate_similarities[2],calculate_similarities[2],calculate_similarities[2],calculate_similarities[2]] 
                
                 elif firstBool== True and secndBool == False:
                     pdbFirstSeqInfo = dictionary_chains.get(pdb_list)
                     if pdbFirstSeqInfo is None: continue
                     list_val = []
                     for mm in range(len(pdb_list_scnd)):
                         if pdb_list_scnd[mm] in dictionary_chains:
                             pdbSecondSeqInfo = dictionary_chains.get(pdb_list_scnd[mm])
                             if pdbSecondSeqInfo is None: continue
                         else: 
                             continue
                         calculate_similarities = find_the_similarites(pdbFirstSeqInfo,pdbSecondSeqInfo)
                         list_val.append(calculate_similarities)
                     processed = divide_similarity_wise(list_val)
                     row = [drugid_list[i],drugid_list[j]] + processed[0]+processed[1]+processed[2]
                 elif firstBool== False and secndBool == True:
                     pdbSecondSeqInfo = dictionary_chains.get(pdb_list_scnd)
                     if pdbSecondSeqInfo is None: continue
                     list_val = []
                     for mm in range(len(pdb_list)):
                         pdbFirstSeqInfo = dictionary_chains.get(pdb_list[mm])
                         if pdbFirstSeqInfo is None: continue
                         calculate_similarities = find_the_similarites(pdbFirstSeqInfo,pdbSecondSeqInfo)
                         list_val.append(calculate_similarities)
                     processed = divide_similarity_wise(list_val)
                     row = [drugid_list[i],drugid_list[j]] + processed[0]+processed[1]+processed[2]
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
                     processed = divide_similarity_wise(list_val)
                     row = [drugid_list[i],drugid_list[j]] + processed[0]+processed[1]+processed[2]
                     print(row)
                     writer.writerow(row)