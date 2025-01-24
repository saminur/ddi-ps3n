# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 23:20:02 2020

@author: Saminur Islam
"""



import sys
import argparse
import codecs
from collections import Counter
from lxml import etree



def run(input):
    """Writes relational database text files for drugs, drug_target, drug_target_action, and targets tables.
    Input: the path of a DrugBank xml file.
    Output: text files that can be used as inputs to SQL tables.
    Method:
        1. Read and parse xml file.
        2. Extract data and save records as key, value pairs.
        3. Write output files."""
    #output file names
    drugs_out = input + '.drugs.txt'
    drug_target_out = input + '.drug_target.txt'
    drug_target_action_out = input + '.drug_target_action.txt'
    targets_out = input + '.targets.txt'
    pathway_out = input + '.pathways.txt'
    drugproducts_out=input + '.drug_products.txt'
    drug_interactions_out=input + '.drug_interactions.txt'
    drug_atc_out=input + '.drug_atc.txt'
    #counter for number of records in each file
    record_counts = Counter()

    #open input file and parse xml
    #get drugbank namspace
    print('Reading and parsing xml file.')
    tree = etree.ElementTree(file=input)
    #a few namespace tricks to make the code more readable
    ns = tree.getroot().nsmap
    ns['db'] = ns[None]
    del ns[None]

    drugs = tree.xpath('db:drug', namespaces=ns)

    #define tables
    #key and value name tuples for each table
    drugs_key_names = ('drug_id',)
    drugs_value_names = ['drug_name', 'restrictions', 'drug_description', 'pathways', 'general_references', ]
    drug_target_key_names = ('drug_id', 'target_id', )
    drug_target_value_names = ['drug_target_references', ]
    drug_target_action_key_names = ('drug_id', 'target_id', 'action', )
    drug_target_action_value_names = []
    targets_key_names = ('target_id', )
    targets_value_names = ['target_name', ]
    pathway_key_names = ('pathway_id')
    pathway_value_names = ['pathway_name','pathway_category','related_drugs','uniport_ids',]
    drugproduct_key_names = ('id',)
    drugproduct_value_names = ['drug_id','product_name','product_labeller','dosage-form','strength','route','approved','country','source',]
    drug_interactions_key_names = ('id',)
    drug_interactions_value_names = ['drug_bank_id1','drug_bank_id2','drug1_name','drug2_name','description',]
    drug_atc_key_names = ('id',)
    drug_atc_value_names = ['atc_code','level1','name1','level2','name2','level3','name3','level4','name4']
    #dicts for each table
    #in the form {<key tuple>: <value dict>}
    drugs_records = {}
    drug_target_records = {}
    drug_target_action_records = {}
    targets_records = {}
    drug_pathways = {}
    drugproducts = {}
    drug_interactions={}
    drug_atc_codes={}
    #process drug records
    print('Processing records.')
#    counts=1
    countsInteraction=1
    countsProduct=1
    for drug in drugs:
        #initialize dict to save unique records to print at end of run
        #in the form {(<key tuple>): {<value dict>}]}
        drugs_record = {}

        drug_id = [i for i in drug.xpath('db:drugbank-id', namespaces=ns) if i.attrib.get('primary') == 'true'][0].text
        drug_name = drug.xpath('db:name', namespaces=ns)[0].text
        drugs_record['drug_name'] = drug.xpath('db:name', namespaces=ns)[0].text
        drug_description = drug.xpath('db:description', namespaces=ns)[0].text
        #deal with drugs that have no description or stray newlines, linefeeds, and tabs
        drugs_record['drug_description'] = drug_description.strip().replace('\n','').replace('\r', '').replace('\t', ' ') if drug_description else ''
        #restrictions
        drugs_record['restrictions'] = ','.join([g.text for g in drug.xpath('db:groups/db:group', namespaces=ns)])
        #pathways
        drugs_record['pathways'] = ','.join([p.text for p in drug.xpath('db:pathways/db:pathway/db:name', namespaces=ns)])
        
        #general_references
        general_references = drug.xpath('db:general-references', namespaces=ns)[0].text
        drugs_record['general_references'] = general_references.strip().replace('\n','').replace('\r', '').replace('\t', ' ') if general_references else ''
        #add this record to the table
        drugs_records[(drug_id,)] = drugs_record
        record_counts['drugs'] += 1
        
        drugs_atc_record={}

        if len(drug.xpath('db:atc-codes/db:atc-code',namespaces=ns))>0:
            atc=drug.xpath('db:atc-codes/db:atc-code',namespaces=ns)[0]
            drugs_atc_record["atc_code"] = atc.attrib.get('code')
            count=0
            for ii in atc.xpath('db:level', namespaces=ns):
                if count == 0:
                    drugs_atc_record["level1"] = ii.attrib.get('code')
                    drugs_atc_record["name1"] = ii.text
                    count=count+1
                    continue
                if count == 1:
                    drugs_atc_record["level2"] = ii.attrib.get('code')
                    drugs_atc_record["name2"] = ii.text
                    count=count+1
                    continue
                if count == 2:
                    drugs_atc_record["level3"] = ii.attrib.get('code')
                    drugs_atc_record["name3"] = ii.text
                    count=count+1
                    continue
                if count == 3:
                    drugs_atc_record["level4"] = ii.attrib.get('code')
                    drugs_atc_record["name4"] = ii.text
                    count=count+1
                
        
            drug_atc_codes[(drug_id,)] = drugs_atc_record
            record_counts['atc-code'] += 1
        for drugproduct in drug.xpath('db:products/db:product',namespaces=ns):
            drugproductrecord={}
            drugproductrecord['drug_id'] = str(drug_id)
            drugproductrecord['product_name'] = drugproduct.xpath('db:name',namespaces=ns)[0].text
            drugproductrecord['product_labeller'] = drugproduct.xpath('db:labeller',namespaces=ns)[0].text
#            drugproductrecord['ema-product-code'] = drugproduct.xpath('db:ema-product-code',namespaces=ns)[0].text
#            drugproductrecord['ema-ma-number'] = drugproduct.xpath('db:ema-ma-number',namespaces=ns)[0].text
#            drugproductrecord['started-marketing-on'] = drugproduct.xpath('db:started-marketing-on',namespaces=ns)[0].text
#            drugproductrecord['end-marketing-on'] = drugproduct.xpath('db:ended-marketing-on',namespaces=ns)[0].text
#
            drugproductrecord['dosage-form'] = drugproduct.xpath('db:dosage-form',namespaces=ns)[0].text
            drugproductrecord['strength'] = drugproduct.xpath('db:strength',namespaces=ns)[0].text
            drugproductrecord['route'] = drugproduct.xpath('db:route',namespaces=ns)[0].text
            drugproductrecord['approved'] = drugproduct.xpath('db:approved',namespaces=ns)[0].text
            drugproductrecord['country'] = drugproduct.xpath('db:country',namespaces=ns)[0].text
            drugproductrecord['source'] = drugproduct.xpath('db:source',namespaces=ns)[0].text
            drugproducts[(str(countsProduct)),] = drugproductrecord
            countsProduct=countsProduct+1
            record_counts['drug products']+=1
            
        for druginteraction in drug.xpath('db:drug-interactions/db:drug-interaction',namespaces=ns):
            drug_interaction_record={}
            drug_interaction_record['drug_bank_id1'] = str(drug_id)
            drug_interaction_record['drug_bank_id2'] = druginteraction.xpath('db:drugbank-id',namespaces=ns)[0].text
            drug_interaction_record['drug1_name'] = drug_name
            drug_interaction_record['drug2_name'] = druginteraction.xpath('db:name',namespaces=ns)[0].text
            drug_interaction_record['description'] = druginteraction.xpath('db:description',namespaces=ns)[0].text
            drug_interactions[(str(countsInteraction)),] = drug_interaction_record
            countsInteraction=countsInteraction+1
            record_counts['drug interactions']+=1
        
                        
        for pathway in drug.xpath('db:pathways/db:pathway',namespaces=ns):
            pathway_record ={}
            pathway_id=pathway.xpath('db:smpdb-id',namespaces=ns)[0].text
#            pathway_record['pathway_id']=pathway_id
            pathway_record['pathway_name'] = pathway.xpath('db:name',namespaces=ns)[0].text
            pathway_record['pathway_category'] = pathway.xpath('db:category',namespaces=ns)[0].text
            pathway_record['related_drugs'] = ','.join([p.text for p in pathway.xpath('db:drugs/db:drug/db:drugbank-id', namespaces=ns)])
            pathway_record['uniport_ids'] = ','.join([p.text for p in pathway.xpath('db:enzymes/db:uniprot-id', namespaces=ns)])
            drug_pathways[(pathway_id),] = pathway_record
            record_counts['pathways'] += 1

        #process the targets of each drug
        for target in drug.xpath('db:targets/db:target', namespaces=ns):
            targets_record = {}
            drug_target_record = {}
            target_id = target.xpath('db:id', namespaces=ns)[0].text
            targets_record['target_name'] = target.xpath('db:name', namespaces=ns)[0].text
            targets_records[(target_id,)] = targets_record
            record_counts['targets'] += 1
            #drug_target_references
            drug_target_references = target.xpath('db:references', namespaces=ns)[0].text
            drug_target_record['drug_target_references'] = drug_target_references.strip().replace('\n', '').replace('\r', '').replace('\t',' ') if drug_target_references else ''
            drug_target_records[(drug_id, target_id,)] = drug_target_record
            record_counts['drug_target'] += 1
#
#            #process the actions of the drug on the target
            for action in drug.xpath('db:targets/db:target/db:actions/db:action', namespaces=ns):
                drug_target_action_record = {}
                action = action.text
                drug_target_action_records[(drug_id, target_id, action,)] = drug_target_action_record
                record_counts['drug_target_action'] += 1

    #open output files files, using codecs to handle data such as 'Dihomo-γ-linolenic acid'
    print('Writing output files.')
    with codecs.open(drugs_out, 'w', encoding='utf-8') as drugs_f, \
            codecs.open(drug_target_out, 'w', encoding='utf-8') as drug_target_f, \
            codecs.open(drug_target_action_out, 'w', encoding='utf-8') as drug_target_action_f, \
            codecs.open(targets_out, 'w', encoding='utf-8') as targets_f, \
            codecs.open(pathway_out, 'w', encoding='utf-8') as pathways_f, \
            codecs.open(drug_interactions_out, 'w', encoding='utf-8') as drug_interactions_f, \
            codecs.open(drugproducts_out, 'w', encoding='utf-8') as drugproducts_f, \
            codecs.open(drug_atc_out, 'w', encoding='utf-8') as drugatccode_f:

        #write output files from the dicts by iteration over each table, output file pair
        for table, file, key_names, value_names in (
                (drugs_records, drugs_f, drugs_key_names, drugs_value_names),
                (targets_records, targets_f, targets_key_names, targets_value_names),
                (drug_target_records, drug_target_f, drug_target_key_names, drug_target_value_names),
                (drug_target_action_records, drug_target_action_f, drug_target_action_key_names, drug_target_action_value_names),
                (drug_pathways, pathways_f, pathway_key_names, pathway_value_names),
                (drug_interactions, drug_interactions_f, drug_interactions_key_names, drug_interactions_value_names),
                (drugproducts, drugproducts_f, drugproduct_key_names, drugproduct_value_names),
                (drug_atc_codes, drugatccode_f, drug_atc_key_names, drug_atc_value_names),
        ):
            #write output file column header
            k = '\t'.join(key_names)
            tab = '\t' if value_names else ''
            v = '\t'.join(value_names)
            file.write(u'{keys}{tab}{values}\n'.format(keys=k,  tab=tab,  values=v))

            #write a key, value record
            for key, value in table.items():
                k = '\t'.join(key)
                tab = '\t' if value_names else ''
                cc=[]
                for n in value_names:
                    if value[n] == None:
                        cc.append('')
                    else:
                        cc.append(value[n])
                    
                v = '\t'.join(cc)
                file.write(u'{keys}{tab}{values}\n'.format(keys=k, tab=tab, values=v))

    print('Done:\n{}'.format('\n'.join(['\t{}: {} records'.format(c, record_counts[c]) for c in sorted(record_counts.keys())])))


path = "full database.xml"
def main():
    run(path)


if __name__ == "__main__": sys.exit(main())