import os
import pandas as pd
from collections import defaultdict

class ProteinKinaseParser:
    """
    Simple class to parse the input/prot_kinase.tsv file and the tdark_kinase.tsv file
    """
    def __init__(self) -> None:
        d = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) 
        self._prot_kinase_tsv_path = os.path.join(d, 'input', 'prot_kinase.tsv')
        if not os.path.exists(self._prot_kinase_tsv_path):
            raise FileNotFoundError("Could not find file at %s" % self._prot_kinase_tsv_path)
        self._dark_kinase_tsv_path = os.path.join(d, 'input', 'tdark_kinase.tsv')
        if not os.path.exists(self._dark_kinase_tsv_path):
            raise FileNotFoundError("Could not find file at %s" % self._dark_kinase_tsv_path)

    def get_symbol_to_id_map(self):
        """
        
        returns a dictionary like this: {'NCBIGene:2870': 'GRK6', 'NCBIGene:140609': 'NEK7', ... }
        """
        print("[INFO] Reading protein kinase information from %s" % self._prot_kinase_tsv_path)
        kinase_data = pd.read_csv(self._prot_kinase_tsv_path, sep="\t", header=None)
        symbol_to_id_map = defaultdict(str)
        for i in range(kinase_data.shape[0]):
            gene_symbol = kinase_data.iloc[i][0]
            ncbi_id = kinase_data.iloc[i][2]
            ncbigene_id = "ncbigene%d" % ncbi_id
            symbol_to_id_map[gene_symbol] =  ncbigene_id
        return symbol_to_id_map

    def get_dark_kinase_map(self):
        print("[INFO] Reading dark kinase information from %s" % self._dark_kinase_tsv_path)
        dark_kinase_data = pd.read_csv(self._dark_kinase_tsv_path, sep="\t")
        symbol_to_id_map = defaultdict(str)
        print(dark_kinase_data.head())
        for i in range(dark_kinase_data.shape[0]):
            gene_symbol = dark_kinase_data.iloc[i][0]
            ncbi_id = dark_kinase_data.iloc[i][1]
            ncbigene_id = "ncbigene%d" % ncbi_id
            symbol_to_id_map[gene_symbol] =  ncbigene_id
        return symbol_to_id_map

    def get_dark_kinase_df(self) -> pd.DataFrame:
        dark_kinase_map = self.get_dark_kinase_map()
        entries = []
        for k,v in dark_kinase_map.items():
            d = {'symbol': k, 'gene_id': v}
            entries.append(d)
        df =  pd.DataFrame(entries)
        cols = ['symbol', 'gene_id']
        return df[cols]
