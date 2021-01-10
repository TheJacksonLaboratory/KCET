import os
from collections import defaultdict


class PkiToKinaseParser:
    """
    Simple class to parse the input/drug_kinase_links.tsv file
    """
    def __init__(self):
        d = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) 
        drug_kinase_links_tsv_path = os.path.join(d, 'input', 'drug_kinase_links.tsv')
        if not os.path.exists(drug_kinase_links_tsv_path):
            raise FileNotFoundError("Could not find file at %s" % drug_kinase_links_tsv_path)
        self._kinase_links_path = drug_kinase_links_tsv_path

    def get_pki_to_kinase_list_dict(self):
        """
        Create a dictionary with the data from drug_kinase_links.tsv
        key -- a protein kinase inhibitor such as abemaciclib	
        value -- list of kinases inhibited by the PKI, e.g., [CDK4,CDK6]
        """
        pki_to_kinase = defaultdict(list)
        with open(self._kinase_links_path) as f:
            for line in f:
                fields = line.rstrip().split('\t')
                if len(fields) != 3:
                    raise ValueError("Bad line in %s (%s)" % (path, line))
                pki = fields[0]
                pk = fields[1] # kinase that is inhibited by the PKI in fields[0]
                pki_to_kinase[pki].append(pk)
        return pki_to_kinase
