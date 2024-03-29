from .drugcentral_pk_pki_parser import DrugCentralPkPkiParser
import os
import pandas as pd
import numpy as np
from typing import List, Dict, Set
from collections import defaultdict
import logging

logging.basicConfig(filename='kcet.log', level=logging.INFO)


class KcetParser:
    """
    Simple class to parse various files needed throughout KCET including input/prot_kinase.tsv file and 
    the tdark_kinase.tsv file
    """

    def __init__(self) -> None:
        """
        Initilized file paths from the ``input`` subfolder of the directory 
        Input the various files into data structures
        """
        # Check that we can find all of the files we need before we start.
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
        # Ingest data
        self._symbol_to_id_map = self._ingest_symbol_to_id_map()
        logging.info("ingested symbol_to_id_map with %d entries such as {'NCBIGene:2870': 'GRK6'}" % len(
            self._symbol_to_id_map))
        # Get reverse map
        self._id_to_symbol_map = {v: k for k, v in self._symbol_to_id_map.items()}
        self._mesh_list = self._ingest_mesh_id_list()
        logging.info("Ingested mesh_id list with %d entries such as 'meshd000008' and 'meshd000069293', " % len(
            self._mesh_list))
        self._meshid2disease_map = self._ingest_mesh_to_disease_map()
        logging.info("Ingested _meshid2disease_map with %d entries" % len(self._meshid2disease_map))
        self._sym2tdl = self._ingest_symbol_to_tdl_map()
        logging.info("Ingested meshid2disease_map with %d entries" % len(self._sym2tdl))
        #self._pki_to_kinase = self._ingest_pki_to_kinase_list_dict()
        self._drug_central = DrugCentralPkPkiParser()

    def _ingest_symbol_to_id_map(self) -> Dict:
        """
        returns a dictionary like this: {'NCBIGene:2870': 'GRK6', 'NCBIGene:140609': 'NEK7', ... }
        """
        logging.info("Reading protein kinase information from %s" % self._prot_kinase_tsv_path)
        kinase_data = pd.read_csv(self._prot_kinase_tsv_path, sep="\t", header=None)
        symbol_to_id_map = defaultdict(str)
        for i in range(kinase_data.shape[0]):
            gene_symbol = kinase_data.iloc[i][0]
            ncbi_id = kinase_data.iloc[i][2]
            ncbigene_id = "ncbigene%d" % ncbi_id
            symbol_to_id_map[gene_symbol] = ncbigene_id
        return symbol_to_id_map

    def get_symbol_to_id_map(self) -> Dict:
        return self._symbol_to_id_map

    def get_id_to_symbol_map(self) -> Dict:
        return self._id_to_symbol_map

    def get_mesh_id_list(self) -> List:
        return self._mesh_list

    def _ingest_mesh_id_list(self) -> List:
        """
        Return a list of the MeSH ids that correspond to cancers from the 
        inputs/neoplasms_labels,tsv file
        """
        mesh_list = []
        with open(self._neoplasms_labels_tsv_path) as f:
            for line in f:
                fields = line.rstrip().split('\t')
                if len(fields) < 2 or len(fields) > 3:
                    # At a minimum, there is MeSH id and label. Most lines have a third field (synonyms)
                    raise ValueError("Bad line in neoplasms_labels,tsv file: %s" % line)
                mesh = fields[0]
                mesh_id_first_letter = mesh[0].lower()
                mesh_id = "mesh" + mesh_id_first_letter + mesh[1:]
                mesh_list.append(mesh_id)
        return mesh_list

    def get_mesh_to_disease_map(self) -> Dict:
        return self._meshid2disease_map

    def _ingest_mesh_to_disease_map(self) -> Dict:
        """
        return a map with 'meshd000008': 'Abdominal Neoplasms', 'meshd000069293': 'Plasmablastic Lymphoma', ...
        """
        meshid2disease_map = defaultdict()
        with open(self._neoplasms_labels_tsv_path) as f:
            for line in f:
                fields = line.rstrip().split('\t')
                if len(fields) < 2 or len(fields) > 3:
                    # At a minimum, there is MeSH id and label. Most lines have a third field (synonyms)
                    raise ValueError("Bad line in neoplasms_labels,tsv file: %s" % line)
                mesh = fields[0]
                mesh_id_first_letter = mesh[0].lower()
                mesh_id = "mesh" + mesh_id_first_letter + mesh[1:]
                disease = fields[1]
                meshid2disease_map[mesh_id] = disease
        return meshid2disease_map

    def decode_predictions(self, vectors: pd.DataFrame, probabilities: np.ndarray,
                           deleteEmbeddings: bool = True) -> pd.DataFrame:
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
        vectors.insert(0, "gene_symbol", gene_symbol_list, True)
        vectors.insert(1, "cancer", cancer_list, True)
        vectors.insert(2, "probability", probabilities, True)
        sorted_vectors = vectors.sort_values(by=['probability'], ascending=False)
        if deleteEmbeddings:
            return sorted_vectors[["gene_symbol", "cancer", "probability"]]
        else:
            return sorted_vectors

    def read_target_level_df(self) -> pd.DataFrame:
        """
        Read the target development level file that related proteins to the levels
        Tdark, Tchem, Tbio, Tclin
        """
        predictions = pd.read_csv(self._target_level_tsv_path, sep="\t")
        columns = ['Sym', 'UniProt', 'Description', 'GeneID', 'TDL']
        return predictions[columns]

    def _ingest_symbol_to_tdl_map(self) -> Dict:
        predictions = self.read_target_level_df()
        sym2tdl = defaultdict(str)
        for _, row in predictions.iterrows():
            symbol = row['Sym']
            tdl = row['TDL']
            sym2tdl[symbol] = tdl
        return sym2tdl

    def get_symbol_to_tdl_map(self) -> Dict:
        return self._sym2tdl




