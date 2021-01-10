import os
import pandas as pd
from collections import defaultdict

class ProteinKinaseParser:
    """
    Simple class to parse the input/prot_kinase.tsv file
    """
    def __init__(self) -> None:
        d = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) 
        self._prot_kinase_tsv_path = os.path.join(d, 'input', 'prot_kinase.tsv')
        if not os.path.exists(self._prot_kinase_tsv_path):
            raise FileNotFoundError("Could not find file at %s" % self._prot_kinase_tsv_path)
        print("[INFO] Reading protein kinase information from %s" % self._prot_kinase_tsv_path)

    def get_symbol_to_id_map(self):
        """
        path: path to the 'input/prot_kinase.tsv file
        returns a dictionary like this: {'NCBIGene:2870': 'GRK6', 'NCBIGene:140609': 'NEK7', ... }
        """
        kinase_data = pd.read_csv(self._prot_kinase_tsv_path, sep="\t", header=None)
        symbol_to_id_map = defaultdict(str)
        for i in range(kinase_data.shape[0]):
            gene_symbol = kinase_data.iloc[i][0]
            ncbi_id = kinase_data.iloc[i][2]
            ncbigene_id = "ncbigene%d" % ncbi_id
            symbol_to_id_map[gene_symbol] =  ncbigene_id
        return symbol_to_id_map