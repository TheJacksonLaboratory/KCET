import random
import pandas as pd
import datetime
from pandas.core.frame import DataFrame
from typing import List

from pandas.io.sql import DatabaseError

from .pki_to_kinase_parser import PkiToKinaseParser
from .protein_kinase_parser import ProteinKinaseParser
from .ct_by_phase_parser import CTParserByPhase
from .neoplasm_parser import NeoplasmParser



class Link:
    """
    Simple class to keep track of positive and ne             gative links using a Hash
    """
    def __init__(self, kinase:str, cancer:str) -> None:
        self._kinase = kinase
        self._cancer = cancer
        
    def __hash__(self):
        return hash((self._kinase, self._cancer))

    def __eq__(self, other):
        return (self.__class__ == other.__class__ and self._kinase == other._kinase and self._cancer == other._cancer)

    def to_dict(self):
        return { 'mesh_id': self._cancer, 'gene_id': self._kinase}

    @staticmethod
    def fromDataFrameToLinkSet(df : pd.DataFrame) -> List:
        linkset = set()
        for index, row in df.iterrows():
            m = row['mesh_id']
            g = row['gene.id']
            L = Link(kinase=g, cancer=m)
            linkset.add(L)
        return linkset



class TestTrainingPredictionGenerator:
    """
    Class to generate test, training, and prediction files.
    """
    def __init__(self, clinical_trials: str, year: int = None) -> None: 
        pkiParser = PkiToKinaseParser()
        self._pki_to_kinase_dict = pkiParser.get_pki_to_kinase_list_dict()
        proteinParser = ProteinKinaseParser()
        self._symbol_to_id_map = proteinParser.get_symbol_to_id_map()
        if year is None:
            # If user does not pass year, take the current year
            now = datetime.datetime.now()
            year = int(now.year)
        self._year = year
        parser = CTParserByPhase(clinical_trials=clinical_trials, year=year)
        self._df_allphases = parser.get_all_phases_for_training() # all positive data, phasae 1,2,3,4
        self._df_phase4 = parser.get_phase_4() # all positive data, phase 4 only
        nparser = NeoplasmParser()
        self._mesh_list = nparser.get_mesh_id_list()


    def get_positive_data_set(self) -> pd.DataFrame:
        """
        Positive training set: all links before the year given in the constructor
        """
        return self._df_phase4

    def write_positive_dataset(self, outfilename: str) -> None:
        self._df_phase4.to_csv(outfilename, sep='\t')

    def get_negative_training_dataset(self, factor:int=10) -> pd.DataFrame:
        """
        Get a negative training set.
        We take Random non-links that were not listed in any of phase 1,2,3,4 in the year up 
        to and including self._year
        We return factor times as many negative examples as we have positive examples
        """
        n_pos_examples = len(self._df_phase4)
        n_neg_examples = n_pos_examples * factor
        kinase_list = [ geneid for _, geneid in self._symbol_to_id_map.items()]
        cancer_id_list = self._mesh_list
        positive_links = Link.fromDataFrameToLinkSet(self._df_allphases)
        negative_links = set()
        i = 0  # use i to limit the number of attempts in case there is some problem
        while len(negative_links) < n_neg_examples and i < 1e6:
            i += 1
            random_cancer = random.choice(cancer_id_list)
            random_kinase = random.choice(kinase_list)
            randomLink = Link(kinase=random_kinase, cancer=random_cancer)
            if randomLink in positive_links:
                print("Skipping random link %s since we found it in the positive set")
                continue
            if randomLink in negative_links:
                print("Skipping random link %s since we already added it to the negative set")
                continue
            negative_links.add(randomLink)
        print("[INFO] We generated a negative set with %d examples (the positive set has %d)" % (len(negative_links), len(positive_links)))
        # format as a pandas dataframe
        negative_dict_list = [ l.to_dict() for l in negative_links]
        return pd.DataFrame(negative_dict_list)

    def write_negative_dataset(self, outfilename: str, factor:int = 10) -> None:
        df = self.get_negative_training_dataset(factor=factor)
        df.to_csv(outfilename, sep='\t')
        

    def get_prediction_dataset(self) -> pd.DataFrame:
        """
        This set contains all possible links between untargeted protein kinases
        and the cancers except the positive dataset calculated up to self._year
        """
        positive_links = Link.fromDataFrameToLinkSet(self._df_allphases)
        kinase_list = [ geneid for _, geneid in self._symbol_to_id_map.items()]
        cancer_id_list = self._mesh_list
        prediction_list = []
        for kinase in kinase_list:
            for cancer in cancer_id_list:
                L = Link(cancer=cancer, kinase=kinase)
                if L in positive_links:
                    continue # Do not include positive links in the prediction set
                prediction_list.append(L.to_dict())
        return pd.DataFrame(prediction_list)

    def write_prediction_dataset(self, outfilename: str) -> None:
        df = self.get_prediction_dataset()
        df.to_csv(outfilename, sep='\t')

    def get_year(self):
        return self._year

        


       
