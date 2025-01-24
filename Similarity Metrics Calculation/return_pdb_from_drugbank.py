# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 01:56:06 2020

@author: Saminur Islam
"""


import csv

drugbank_pdb_dir = 'drugTable.csv'
new_pdb_list = []
drugid_list = []
dictionary = dict()

def return_the_desired_list(pdbid_list):
    with open(drugbank_pdb_dir,'r') as drugbankFile:
        count = 0
        readCSV = csv.reader(drugbankFile)
        for row in readCSV:
            if count == 0: 
                count=count+1
                continue
            if len(row) < 6: continue 
            pdb1=row[6]
            if len(row) < 8: 
                pdb2=''
             
            else:
                pdb2=row[8]
            if len(row) < 10: 
                pdb3=''
             
            else:
                pdb3=row[10]    
            DrugId = row[2]  
            drugid_list.append(DrugId)
            if dictionary.get(DrugId):
                listVal=dictionary.get(DrugId)
                if type(listVal)==str:
                    if pdb1 in pdbid_list:
                        plusVal=[listVal]+[pdb1]
                        new_pdb_list.append(pdb1)
                    elif pdb2 in pdbid_list:
                        plusVal=[listVal]+[pdb2]
                        new_pdb_list.append(pdb2)
                    elif pdb3 in pdbid_list:
                        plusVal=[listVal]+[pdb3]
                        new_pdb_list.append(pdb3)
                
                else:
                    if pdb1 in pdbid_list:
                        plusVal=listVal+[pdb1]
                        new_pdb_list.append(pdb1)
                    elif pdb2 in pdbid_list:
                        plusVal=listVal+[pdb2]
                        new_pdb_list.append(pdb2)
                    elif pdb3 in pdbid_list:
                        plusVal=listVal+[pdb3]   
                        new_pdb_list.append(pdb3)
                removeDupl=list(dict.fromkeys(plusVal))
                dictionary[DrugId]=removeDupl
            else:
                if pdb1 in pdbid_list:
                    dictionary[DrugId]=pdb1
                    new_pdb_list.append(pdb1)
                elif pdb2 in pdbid_list:
                    dictionary[DrugId]=pdb2
                    new_pdb_list.append(pdb2)
                elif pdb3 in pdbid_list:
                    dictionary[DrugId]=pdb3
                    new_pdb_list.append(pdb3)
    return new_pdb_list,dictionary,drugid_list