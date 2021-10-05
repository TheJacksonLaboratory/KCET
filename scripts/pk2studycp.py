import argparse
import os
import sys
import csv

sys.path.insert(0, os.path.abspath('..'))

# protein kinase (PK) to study cross product (pk2studycp.py)
# extract a list of all indications tested for protein kinase inhibitors with activity against a specific PK

from kcet import PkPkiFilter, CTParserByPhase, KcetDatasetGenerator

parser = argparse.ArgumentParser(description='Process PKI/PK data')
parser.add_argument('--pk', type=str, default='FLT3')
parser.add_argument('--clinicaltrials', default='/home/peter/data/pubmed2vec/clinical_trials_by_phase.tsv')
parser.add_argument('--outfilename', type=str, default='drug_kinase_links.tsv')
args = parser.parse_args()

the_pk = args.pk
the_ctfile = args.clinicaltrials

print("Extraciong all indications tested for protein kinase inhibitors with activity against {}".format(the_pk))

pkpki = PkPkiFilter()
all_pki = set()

for _, row in pkpki.get_all_pk_pki().iterrows():
    all_pki.add(row['PKI'])

print("Total of {} PKIs".format(len(all_pki)))

disease_indications = set()

with open(the_ctfile) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        disease = row['disease']
        drug = row['drug']
        phase = row['phase']
        if phase == 'Phase 4':
            if drug in all_pki:
                disease_indications.add(disease)

for disease in sorted(list(disease_indications)):
    print(disease)

download_dir = '/home/peter/data/pubmed2vec'
embeddings2020 = os.path.join(download_dir, "embedding_SG_dim100_upto2020.npy")
words2020 = os.path.join(download_dir, "words_SG_upto2020.txt")
ctrials = os.path.join(download_dir,"clinical_trials_by_phase.tsv")

datagen = KcetDatasetGenerator(clinical_trials=ctrials, embeddings=embeddings2020,words=words2020)
X = datagen.get_all_phases_all_pk_pki(target_year=2020)
print("X, ", type(X), "size = ", len(X))

for x in X:
    print(x)
