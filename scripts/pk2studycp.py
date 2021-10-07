import argparse
import os
import sys

sys.path.insert(0, os.path.abspath('..'))

# This script was used to double-check the extraction of PK-cancer links used for testing.

from kcet import DrugCentralPkPkiParser, CTParserByPhase, KcetDatasetGenerator

parser = argparse.ArgumentParser(description='Process PKI/PK data')
parser.add_argument('--pk', type=str, default='FLT3')
parser.add_argument('--clinicaltrials', default='/home/peter/data/pubmed2vec/clinical_trials_by_phase.tsv')
parser.add_argument('--embeddings', default='/home/peter/data/pubmed2vec/embedding_SG_dim100_upto2010.npy')
parser.add_argument('--words', default='/home/peter/data/pubmed2vec/words_SG_upto2010.txt')


parser.add_argument('--outfilename', type=str, default='drug_kinase_links.tsv')
args = parser.parse_args()

the_pk = args.pk
the_ctfile = args.clinicaltrials
embeddings = args.embeddings
words = args.words

print("Extracting all indications tested for protein kinase inhibitors with activity against {}".format(the_pk))

datagen = KcetDatasetGenerator(clinical_trials = the_ctfile,  embeddings=embeddings, words=words)



all_phases = datagen._df_allphases[datagen._df_allphases['year'] <= 2020]
ap_df = all_phases[all_phases['kinase'] == the_pk]
cancer_s = set()
for i, row in ap_df.iterrows():
    cancer_s.add(row['cancer'])
for idx, cancer in enumerate(sorted(list(cancer_s))):
    print("{}) {} ({})".format(idx, cancer, the_pk))
