# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 15:16:07 2020

@author: Saminur Islam
"""


import requests
import time
from pandas import read_excel

my_sheet = 'Sheet' # change it to your sheet name
file_name = 'uniprot.xlsx' # change it to the name of your excel file
df = read_excel(file_name, sheet_name = my_sheet)
uniprot_ids=[]

for i in range(len(df)):
    uid = df.iloc[i][0]
    uniprot_ids.append(uid)
            

uniprot_ids=list(dict.fromkeys(uniprot_ids))


# uniprot_ids = ['P26378', 'O35433', 'Q02910']
url = 'https://www.uniprot.org/uniprot/'

protein_to_pdb = {}
for protein in uniprot_ids:
    params = {
        'format': 'tab',
        'query': 'ID:{}'.format(protein),
        'columns': 'id,database(PDB)'
    }
    contact = ""  # Please set your email address here.
    headers = {'User-Agent': 'Python {}'.format(contact)}
    r = requests.get(url, params=params, headers=headers)

    protein_to_pdb[protein] = str(r.text).splitlines()[-1].split('\t')[-1].split(';')
    protein_to_pdb[protein].pop(-1)
    time.sleep(1)  # be respectful and don't overwhelm the server with requests

test_dictionary = protein_to_pdb
d = { k : v for k,v in protein_to_pdb.items() if v}
saving_directory='uniprot_pdb.csv'
with open(saving_directory, 'w',newline='') as f:
    writer = csv.writer(f)
    for key in d.keys():
        row=[key]+d[key]
        writer.writerow(row)
print(protein_to_pdb)  


import csv
from pandas import read_excel

uniprot_pdb_path = 'uniprot_pdb.csv'
uniprot_dict = dict() 

def function_uniprot_drug_pdb_relation():
    with open(uniprot_pdb_path) as csv_file:
        read_csv = csv.reader(csv_file)
        count = 0
        for row in read_csv:
            uniprot_dict[row[0]] = row[1:]

    my_sheet = 'Sheet1' # change it to your sheet name
    file_name = 'pathways.xlsx' # change it to the name of your excel file
    df = read_excel(file_name, sheet_name = my_sheet)
    dictionary={}
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
    drug_pdb_dict = dict()
    for ii in range(len(drugList)):
        uniprot_idList = dictionary.get(drugList[ii])
        pdb_list =[]
        for jj in range(len(uniprot_idList)):
            if uniprot_idList[jj] in uniprot_dict:
                pdb_values = uniprot_dict.get(uniprot_idList[jj])
                if type(pdb_values) == str:
                    pdb_list.append(pdb_values)
                else:
                    for kk in range(len(pdb_values)):
                        pdb_list.append(pdb_values[kk])
                    
        pdb_list=list(dict.fromkeys(pdb_list))
        drug_pdb_dict[drugList[ii]] = pdb_list
    
    return drug_pdb_dict
    
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
                        if type(uniprot_list) == str:
                            new_list = [uniprot_list] + [uniprot_id]
                        else:
                            new_list = uniprot_list + [uniprot_id]
                    
                        drugbank_uniprot[drug_ids] = new_list
                    else:
                        drugbank_uniprot[drug_ids] = uniprot_id
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
                            drugbank_uniprot[drug_ids[i].replace(" ", "")] = uniprot_id
                        
                    
        
                pdb_ids = str(row[7]).split(";")
                pdb_ids = [x.strip(' ') for x in pdb_ids]
                if type(pdb_ids) == str:
                    uniprot_pdb_map[uniprot_id] = [pdb_ids]
                else:
                    uniprot_pdb_map[uniprot_id] = pdb_ids

    drugbank_pdb_map = dict()
    for drugbankId in drugbank_uniprot: 
        uniprotIds = drugbank_uniprot.get(drugbankId)
        if type(uniprotIds) == str:
           if uniprotIds in uniprot_pdb_map:
               pdbIds = uniprot_pdb_map.get(uniprotIds)
               for m in range(len(pdbIds)):
                   if pdbIds[m] == '': 
                       continue
                   if drugbank_pdb_map.get(drugbankId):
                       pdb_list = drugbank_pdb_map.get(drugbankId)
                       if type(pdb_list) == str:
                           new_list = [pdb_list] + [pdbIds[m]]
                       else:
                           new_list = pdb_list + [pdbIds[m]]
                       drugbank_pdb_map[drugbankId] = new_list
                   else:
                       drugbank_pdb_map[drugbankId] = pdbIds[m]
                    
                    
        else:
            for k in range(len(uniprotIds)):
                if uniprotIds[k] in uniprot_pdb_map:
                    pdbIds = uniprot_pdb_map.get(uniprotIds[k])
                    for m in range(len(pdbIds)):
                        if pdbIds[m] == '': 
                            continue
                        if drugbank_pdb_map.get(drugbankId):
                            pdb_list = drugbank_pdb_map.get(drugbankId)
                            if type(pdb_list) == str:
                                new_list = [pdb_list] + [pdbIds[m]]
                            else:
                                new_list = pdb_list + [pdbIds[m]]
                    
                            drugbank_pdb_map[drugbankId] = new_list
                        else:
                            drugbank_pdb_map[drugbankId] = pdbIds[m]
    return drugbank_pdb_map
            
drugbank_pdb_map = pdbid_drugid_relation()

# for ii in drugbank_pdb_map:
#     pdbids = drugbank_pdb_map.get(ii)
#     if type(pdbids) != str:
#         pdbids = list(dict.fromkeys(pdbids))
#         drugbank_pdb_map[ii] = pdbids
                    
                    