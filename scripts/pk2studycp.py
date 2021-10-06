import argparse
import os
import sys
import csv

sys.path.insert(0, os.path.abspath('..'))

# protein kinase (PK) to study cross product (pk2studycp.py)
# extract a list of all indications tested for protein kinase inhibitors with activity against a specific PK

from kcet import DrugCentralPkPkiParser, CTParserByPhase, KcetDatasetGenerator

parser = argparse.ArgumentParser(description='Process PKI/PK data')
parser.add_argument('--pk', type=str, default='FLT3')
parser.add_argument('--clinicaltrials', default='/home/peter/data/pubmed2vec/clinical_trials_by_phase.tsv')
parser.add_argument('--outfilename', type=str, default='drug_kinase_links.tsv')
args = parser.parse_args()

the_pk = args.pk
the_ctfile = args.clinicaltrials

print("Extracting all indications tested for protein kinase inhibitors with activity against {}".format(the_pk))

parser = CTParserByPhase(clinical_trials=the_ctfile)
df = parser.get_all_phases() # _drug_kinase_links has entries like {abemaciclib:[CDK4,CDK6]}()

print("Total of {} PKIs".format(len(df)))
print(df.head())

df_for_kinase = df[df['kinase']==the_pk]
cancers = df_for_kinase['cancer'].unique()
disease_indications = set()


for idx, cancer in enumerate(sorted(cancers)):
    print("{}) {}".format(idx, cancer))
