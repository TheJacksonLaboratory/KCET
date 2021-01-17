import os
import pandas as pd
import numpy as np
from collections import defaultdict

class KcetParser:
    """
    Simple class to parse various files needed throughout KCET including input/prot_kinase.tsv file and the tdark_kinase.tsv file
    """
    def __init__(self) -> None:
        d = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) 
        self._prot_kinase_tsv_path = os.path.join(d, 'input', 'prot_kinase.tsv')
        if not os.path.exists(self._prot_kinase_tsv_path):
            raise FileNotFoundError("Could not find file at %s" % self._prot_kinase_tsv_path)
        self._neoplasms_labels_tsv_path = os.path.join(d, 'input', 'neoplasms_labels.tsv')
        if not os.path.exists(self._neoplasms_labels_tsv_path):
            raise FileNotFoundError("Could not find file at %s" % self._neoplasms_labels_tsv_path)
        self._target_level_tsv_path = os.path.join(d, 'input', 'target_develop_levels.tsv.csv')
        if not os.path.exists(self._target_level_tsv_path):
            raise FileNotFoundError("Could not find file at %s" % self._target_level_tsv_path)

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

    def get_id_to_symbol_map(self):
        kinase_data = pd.read_csv(self._prot_kinase_tsv_path, sep="\t", header=None)
        id_to_symbol_map = defaultdict(str)
        for i in range(kinase_data.shape[0]):
            gene_symbol = kinase_data.iloc[i][0]
            ncbi_id = kinase_data.iloc[i][2]
            ncbigene_id = "ncbigene%d" % ncbi_id
            id_to_symbol_map[ncbigene_id] =  gene_symbol
        return id_to_symbol_map

    def get_mesh_id_list(self):
        """
        Return a list of the MeSH ids that correspond to cancers from the 
        inputs/neoplasms_labels,tsv file
        """
        mesh_list = []
        with open(self._neoplasms_labels_tsv_path) as f:
            for line in f:
                fields = line.rstrip().split('\t')
                if len(fields) < 2 or len(fields)> 3:
                    #At a minimum, there is MeSH id and label. Most lines have a third field (synonyms)
                    raise ValueError("Bad line in neoplasms_labels,tsv file: %s"  % line)
                mesh = fields[0]
                mesh_id_first_letter = mesh[0].lower()
                mesh_id = "mesh" + mesh_id_first_letter + mesh[1:]
                mesh_list.append(mesh_id)
        return mesh_list

    def get_mesh_to_disease_map(self):
        """
        return a map with 'meshd000008': 'Abdominal Neoplasms', 'meshd000069293': 'Plasmablastic Lymphoma', ...
        """
        meshid2disease_map = defaultdict()
        with open(self._neoplasms_labels_tsv_path) as f:
            for line in f:
                fields = line.rstrip().split('\t')
                if len(fields) < 2 or len(fields)> 3:
                    #At a minimum, there is MeSH id and label. Most lines have a third field (synonyms)
                    raise ValueError("Bad line in neoplasms_labels,tsv file: %s"  % line)
                mesh = fields[0]
                mesh_id_first_letter = mesh[0].lower()
                mesh_id = "mesh" + mesh_id_first_letter + mesh[1:]
                disease = fields[1]
                meshid2disease_map[mesh_id] = disease
        return meshid2disease_map

    def decode_predictions(self, vectors: pd.DataFrame, probabilities: np.ndarray, deleteEmbeddings: bool=True) -> pd.DataFrame:
        """
        We create difference vectors to have labels like this -- ncbigene5599-meshd000074723
        that is, a gene (protein kinase) and a disease
        For output, it is nice to have the labels for the genes and diseases, so we add two columns
        We also add the probability
        ncbigene94-meshd006528 	ACVRL1 	Carcinoma, Hepatocellular 	1.0 	
        """
        if len(vectors) != len(probabilities):
            raise ValueError("The length of the vectors dataframe and the probabilities do not match!")
        meshid2disease_map = self.get_mesh_to_disease_map()
        ncbigene2symbol_map = self.get_id_to_symbol_map()
        gene_symbol_list = []
        cancer_list = []
        for vec in vectors.index:
            fields = vec.split("-")
            ncbi_gene = fields[0]
            mesh_cancer = fields[1]
            gene_symbol = ncbigene2symbol_map[ncbi_gene]
            gene_symbol_list.append(gene_symbol)
            cancer = meshid2disease_map[mesh_cancer]
            cancer_list.append(cancer)
        vectors.insert(0,"gene_symbol", gene_symbol_list, True)
        vectors.insert(1,"cancer", cancer_list, True)
        vectors.insert(2,"probability",probabilities, True)
        sorted_vectors = vectors.sort_values(by=['probability'],ascending=False)
        if deleteEmbeddings:
            return sorted_vectors[["gene_symbol","cancer","probability"]]
        else:
            return sorted_vectors

    def read_target_level_df(self):
        """
        Read the target development level file that related proteins to the levels
        Tdark, Tchem, Tbio, Tclin
        """
        predictions = pd.read_csv(self._target_level_tsv_path,sep="\t")
        columns = ['Sym', 'UniProt', 'Description', 'GeneID', 'TDL']
        return predictions[columns]

    def get_symbol_to_tdl_map(self):
        predictions = self.read_target_level_df()
        sym2tdl = defaultdict(str)
        for _, row in predictions.iterrows():
            symbol = row['Sym']
            tdl = row['TDL']
            sym2tdl[symbol] = tdl
        return sym2tdl
