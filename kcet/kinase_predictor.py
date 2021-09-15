
from numpy.lib.function_base import median
import pandas as pd 
import numpy as np
import os
from collections import defaultdict



class KinasePredictor:

    def __init__(self, embeddings, words) -> None:
        """
        The constructor ingests the prot_kinase.tsv and neoplasm_labels.tsv files from the input directory
        If also ingests the embeddings/words and places them into a pandas Dataframe
        """
        embedding = np.load(embeddings, mmap_mode=None, allow_pickle=False, fix_imports=True, encoding='ASCII')
        word_list = []
        with open(words) as f:
            for line in f:
                word = line[2:-3]
                word_list.append(word)
        self._embeddings_df = pd.DataFrame(data=embedding,index = word_list)
        print("")
        print("[INFO] We ingested %d labeled word vectors from %s and %s" % (len(self._embeddings_df), embeddings, words))
        path = os.path.dirname(os.path.abspath(__file__)) 
        parent_dir = os.path.abspath(os.path.join(path, os.pardir))
        prot_kinase_tsv = os.path.join(parent_dir, "input/prot_kinase.tsv")
        if not os.path.isfile(prot_kinase_tsv):
            raise FileNotFoundError("Could not find prot_kinase.tsv file at %s" % prot_kinase_tsv)
        kinase_gene_id = pd.read_csv(prot_kinase_tsv, sep = "\t", header = None)
        self._ncbigene2symbol_map = defaultdict()
        for i in kinase_gene_id.index:
            gene_symbol = kinase_gene_id.iloc[i][0]
            ncbigene = kinase_gene_id.iloc[i][2]
            ncbigene_id = "ncbigene" + str(ncbigene)
            self._ncbigene2symbol_map[ncbigene_id] = gene_symbol
        print("[INFO] We ingested %d symbol/NCBI gene id mappings from %s" % (len(self._ncbigene2symbol_map), prot_kinase_tsv))
        neoplasm_labels_tsv = os.path.join(parent_dir, "input/neoplasms_labels.tsv")
        disease_mesh = pd.read_csv(neoplasm_labels_tsv,  sep= "\t", header = None)
        self._meshid2disease_map = defaultdict()
        for i in disease_mesh.index:
            mesh = disease_mesh.iloc[i][0]
            mesh_first_letter = mesh[0].lower()
            mesh_id = "mesh" + mesh_first_letter + mesh[1:]
            disease = disease_mesh.iloc[i][1]
            self._meshid2disease_map[mesh_id] = disease
        print("[INFO] We ingested %d meshId/disease mappings from %s" % (len(self._meshid2disease_map), neoplasm_labels_tsv))
        print("")


    def get_words(self):
        return self._embeddings_df.index

    def get_embeddings(self) -> pd.DataFrame :
        """
        return a Pandas dataframe whose index is the words, and whose columns are the dimensions of the embeddings
        """
        return self._embeddings_df

    def get_disease_kinase_difference_vectors(self, examples: pd.DataFrame) -> pd.DataFrame:
        """
        The input is a dataframe with protein kinases (NCBI gene ids) and cancers (MeSH id)
        This function finds the embedded vectors for the genes and cancers, it substracts the
        cancer vector from the kinase vector, and it returns a data frame with these vectors.
        This method assumees that the input dataframe contains columns called gene_id and mesh_id and will
        fail if this is not the case
        """
        unidentified_genes = set()
        unidentified_cancers = set()
        if not "gene_id" in examples.columns:
            raise ValueError("Input dataframe must contain a column called gene_id")
        if not "mesh_id" in examples.columns:
            raise ValueError("Input dataframe must contain a column called mesh_id")
        # if len(examples.columns) != 2:
        #    raise ValueError("Input dataframe must have exactly two columns")
        df = pd.DataFrame(columns = self._embeddings_df.columns)
        total = len(examples.index)
        if total==0:
            raise ValueError("Attempt to get difference vectors from empty data frame")
        i = 0
        for _, row in examples.iterrows(): 
            ncbigene_id = row["gene_id"]
            mesh_id = row["mesh_id"]
            ncbigene_id_embedding = None
            mesh_id_embedding = None
            if ncbigene_id in self._embeddings_df.index:
                ncbigene_id_embedding = self._embeddings_df.loc[ncbigene_id]
            else:
                unidentified_genes.add(ncbigene_id)
            if mesh_id in self._embeddings_df.index:
                mesh_id_embedding = self._embeddings_df.loc[mesh_id]
            else:
                unidentified_cancers.add(mesh_id)
            if ncbigene_id_embedding is not None and mesh_id_embedding is not None:     
                diff_kinase_mesh = np.subtract(ncbigene_id_embedding, mesh_id_embedding)
                #diff_kinase_mesh_list_pos_train.append(diff_kinase_mesh)
                #diff_index_pos_train.append(ncbigene_id + "," + mesh_id)
                label = "%s-%s" % (ncbigene_id, mesh_id)
                df.loc[label] = diff_kinase_mesh
            i += 1
            if i % 10000 == 0 and i > 0:
                print("[INFO] Created %d/%d (%.1f%%) difference vectors" % (i, total, 100.0*i/total))
        print("[INFO] Extracted %s kinase-cancer difference vectors" % len(df))
        print("[INFO]\tInitial data: %d examples" % len(examples))
        print("[INFO]\tCould not identify %d gene ids" % len(unidentified_genes))
        print("[INFO]\tCould not identify %d MeSH ids" % len(unidentified_cancers))
        return df
