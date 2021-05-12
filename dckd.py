
from kcet import PkPkiFilter

## Temporary file to develop PkPkiFilter


pkpki = PkPkiFilter()
outfilename = 'input/drug_kinase_links.tsv'
for x in range(21):
    valid_pk_pki = pkpki.get_valid_pk_pki(n_pki_limit=x)
    print("Threshold: %d interactions; number of PKI-PK links: %d" %(x, len(valid_pk_pki)))



