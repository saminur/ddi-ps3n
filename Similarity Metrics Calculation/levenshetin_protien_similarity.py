# -*- coding: utf-8 -*-
"""
Created on Wed May 20 15:09:37 2020

@author: Saminur Islam
"""

import textdistance
import math
import numpy as np


columns = ['drug id One','drug id two','Levenshtein distance min','Levenshtein distancee max','Levenshtein distance avg', 'Levenshtein distance weighted avg']

similarity_directory = 'drug_protein_levenshtein_similarity.csv'



def calculate_similarity_from_protein_target(firstListProtein,secndListProtein):
    max_simi = 0
    min_simi = 1
    avg = 0
    exp_list = []
    exponetial_avg = 0
    list_val = []
    for ii in range(len(firstListProtein)):
        pro_tar1= firstListProtein[ii]
        for jj in range(len(secndListProtein)):
            pro_tar2 = secndListProtein[jj]
            levn = textdistance.levenshtein.normalized_similarity(pro_tar1, pro_tar2)
            list_val.append(levn)
            if levn > max_simi:
                max_simi = levn
            if levn < min_simi:
                min_simi = levn
            
            avg = avg + levn
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
    avg = avg/ (len(firstListProtein) * len(secndListProtein))
    weight_list  = np.divide(exp_list, exponetial_avg)
    final_weighted_avg = 0
    for i in range(len(list_val)):
        final_weighted_avg = final_weighted_avg + list_val[i] * weight_list[i]
    
    calulted_value_list = [min_simi,max_simi,avg,final_weighted_avg/len(list_val)]
    return calulted_value_list

with open(similarity_directory, 'w',newline='') as f:
    writer = csv.writer(f)
    writer.writerow(columns)
    for i in range(len(drugList)):
        firstListProtein = dictionary.get(drugList[i])
        for j in range(i+1,len(drugList)):
            secndListProtein = dictionary.get(drugList[j])
            similarity_row = calculate_similarity_from_protein_target(firstListProtein,secndListProtein)
            final_row = [drugList[i],drugList[j]]+similarity_row
            writer.writerow(final_row)
            print(final_row)
            
'''levenshtein distance similarity ended'''    