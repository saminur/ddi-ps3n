# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 23:23:55 2023

@author: samin
"""

from pandas import read_excel
import numpy as np
import pandas as pd
import csv 

path_found_ddi = 'H:\\Research Work\\TrinetXDataSet\\New Results for Jounral\\new_found_ddi_ndd.csv'
path_drugs = 'I:\\Study\\Research WVU\\Drug Bank\\drug_products.txt' 

path_similarity = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\NDD Dataset\\ds1_pdb_ndd\\merged final\\merged_protein_pdb_ds1.csv'
path_simi_cosine_pdb = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\NDD Dataset\\\ds1_pdb_ndd\\pdb matrices\\pdb_cosine_avg.csv'
path_simi_cosine_protein = 'I:\\Study\\Research WVU\\PDB\\drug_similarity\\NDD Dataset\\\ds1_protein_ndd\\protein_cosine_avg.csv'
df_drug_info = pd.read_csv(path_drugs, sep="\t",header=0).drop("id", axis=1)
df_path_found_ddi = pd.read_csv(path_found_ddi)
df_path_smilarity = pd.read_csv(path_similarity).set_index('Unnamed: 0')

df_protein_cosine = pd.read_csv(path_similarity).set_index('Unnamed: 0')
df_pdb_cosine = pd.read_csv(path_simi_cosine_pdb).set_index('Unnamed: 0')
dict_drug_band = {}
for index,row in df_drug_info.iterrows():
    # print(row['drug_id'])
    if row['drug_id'] not in dict_drug_band:
        dict_drug_band[row['drug_id']] = row['brand_name']
    # print(df_drug_info[i]['drug id'])
    # break;
    
dict_similarity = {}

path = 'H:\\Research Work\\TrinetXDataSet\\New Results for Jounral\\network_ndd.csv'
columns = ['Drug ID1', 'Drug Name1', 'Drug ID2', 'Drug Name2', 'Interaction','fushion similarity']
df_path_smilarity['DB00115']['DB00136']

with open(path, 'w',newline='') as f:
    writer = csv.writer(f)
    writer.writerow(columns)
    for index,row in df_path_found_ddi.iterrows():
        DrugName1 = dict_drug_band.get(row['drug1'])
        DrugName2 = dict_drug_band.get(row['Drug2'])
        similarity = df_path_smilarity[row['drug1']][row['Drug2']]
        write_row = [row['drug1'],DrugName1,row['Drug2'],DrugName2,row['Interactions'],similarity]
        writer.writerow(write_row)
        

df_network_ndd = pd.read_csv(path)       
