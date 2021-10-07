import argparse
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# extract list of pk pki links
# view PK-PKI links according to PK-per-PKI (n_pk)
# See manuscript for explanation of n_pk

from kcet import DrugCentralPkPkiParser

parser = argparse.ArgumentParser(description='Process PKI/PK data')
parser.add_argument('--max_multiplicity', type=int, default=5)
parser.add_argument('--outfilename',  type=str, default='drug_kinase_links.tsv')
args = parser.parse_args()

n_pk_pki = args.max_multiplicity
outfile_name = "drug_kinase_links_n_pk_{}.tsv".format(n_pk_pki)
fh = open(outfile_name, 'wt')
pkpki = DrugCentralPkPkiParser()
valid_pk_pki = pkpki.get_pk_pki_with_threshold(n_pki_limit=n_pk_pki)
for pki, kinase_list in valid_pk_pki.items():
    for kinase in kinase_list:
        fh.write("{}\t{}\n".format(pki, kinase))
fh.close()