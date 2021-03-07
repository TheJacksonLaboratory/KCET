import random
import pandas as pd
import datetime
from typing import List, Tuple

from .kcet_parser import KcetParser
from .ct_by_phase_parser import CTParserByPhase



class Link:
    """
    Simple class to keep track of positive and negative links using a Hash
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
            g = row['gene_id']
            L = Link(kinase=g, cancer=m)
            linkset.add(L)
        return linkset

class KcetDatasetGenerator:
    """
    Class to generate test, training, and prediction files.
    """
    def __init__(self, clinical_trials: str) -> None: 
        pkiParser = KcetParser()
        self._pki_to_kinase_dict = pkiParser.get_pki_to_kinase_list_dict()
        kcetParser = KcetParser()
        self._symbol_to_id_map = kcetParser.get_symbol_to_id_map()
        parser = CTParserByPhase(clinical_trials=clinical_trials)
        self._df_allphases = parser.get_all_phases() # all positive data, phase 1,2,3,4
        self._df_phase4 = parser.get_phase_4() # all positive data, phase 4 only
        self._mesh_list = kcetParser.get_mesh_id_list()

    def get_latest_data(self):
        now = datetime.datetime.now() # default to current year
        year = now.year
        return year

    def get_data_for_target_year(self, target_year:int, factor:int=10) -> Tuple[pd.DataFrame]:
        """
        Get positive and negative data for the target year indicated by the argument.
        """
        positive_df = self._get_positive_data_set(year=target_year)
        negative_df = self._get_negative_training_dataset(positive_df=positive_df, factor=factor)
        now = datetime.datetime.now() # default to current year
        currentyear = now.year
        if target_year < currentyear: #historical prediction
            positive_validation_df = self._get_positive_validation_data_set(year=target_year)
            negative_validation_df = self._get_negative_validation_data_set(negative_df=negative_df)
            return positive_df, negative_df, positive_validation_df, negative_validation_df
        else: #novel prediction
            prediction_df = self.get_prediction_dataset()
            return positive_df, negative_df, prediction_df

    def get_data_for_target_year_phase_4(self, target_year: int, factor: int = 10) -> Tuple[pd.DataFrame]:
        """
        Get positive and negative data for the target year indicated by the argument.
        """
        positive_df = self._get_positive_data_set(year=target_year)
        negative_df = self._get_negative_training_dataset(positive_df=positive_df, factor=factor)
        now = datetime.datetime.now()  # default to current year
        currentyear = now.year
        if target_year < currentyear:  # historical prediction
            positive_validation_df = self._get_positive_validation_data_set_phase_4(year=target_year)#get phase 4 for the years after the target year
            negative_validation_df = self._get_negative_validation_data_set(negative_df=negative_df)
            return positive_df, negative_df, positive_validation_df, negative_validation_df
        else:  # novel prediction
            prediction_df = self.get_prediction_dataset()
            return positive_df, negative_df, prediction_df

    def _get_positive_data_set(self, year:int) -> pd.DataFrame:
        """
        Positive training set: all links before the year given in the constructor
        """
        if not isinstance(year, int):
            raise ValueError("year must be an integer")

        not_later_than_target_year = self._df_phase4['year'] <= year
        return self._df_phase4[not_later_than_target_year]

    def _get_positive_validation_data_set(self, year:int) -> pd.DataFrame:
        """
        Get all of the positive examples from after the target year-- used for validation
        in historical validation experiments
        """
        print("GPVDS year=" + str(year))
        later_than_target_year = self._df_allphases['year'] > year
        return self._df_allphases[later_than_target_year]

    def _get_positive_validation_data_set_phase_4(self, year: int) -> pd.DataFrame:
        """
        Get all of the positive links (of phase 4) from after the target year-- used for validation
        in historical validation experiments
        """
        print("GPVDS year=" + str(year))
        phase_4_later_than_target_year = self._df_phase4['year'] > year
        return self._df_phase4[phase_4_later_than_target_year]
    
    def _get_negative_validation_data_set(self, negative_df: pd.DataFrame) -> pd.DataFrame:
        """
        Get negative examples from after the target year-- used for validation
        in historical validation experiments. Note that we take examples that are negative
        from the perspective of the current time -- we are taking factor-times more negative
        examples than positive examples, and this function chooses a set that is distinct
        from the set of examples use prior to the target year (the same size as the
        negative_df that is passed to us).
        """
        n_neg_examples = len(negative_df)
        kinase_list = [ geneid for _, geneid in self._symbol_to_id_map.items()]
        cancer_id_list = self._mesh_list
        positive_links = Link.fromDataFrameToLinkSet(self._df_allphases)
        pretarget_negative_links = Link.fromDataFrameToLinkSet(negative_df)
        negative_links = set()
        i = 0  # use i to limit the number of attempts in case there is some problem
        while len(negative_links) < n_neg_examples and i < 1e6:
            i += 1
            random_cancer = random.choice(cancer_id_list)
            random_kinase = random.choice(kinase_list)
            randomLink = Link(kinase=random_kinase, cancer=random_cancer)
            if randomLink in positive_links:
                print("Skipping random link (%s,%s) since we found it in the positive set" % (random_kinase,random_cancer))
                continue
            if randomLink in pretarget_negative_links:
                print("Skipping random link (%s,%s) since we already added it to the negative set" % (random_kinase,random_cancer))
                continue
            negative_links.add(randomLink)
        print("[INFO] We generated a negative set with %d examples (the positive set has %d)" % (len(negative_links), len(positive_links)))
        # format as a pandas dataframe
        negative_dict_list = [ l.to_dict() for l in negative_links]
        return pd.DataFrame(negative_dict_list)

    def _get_negative_training_dataset(self, positive_df: pd.DataFrame, factor:int=10) -> pd.DataFrame:
        """
        Get a negative training set.
        We take Random non-links that were not listed in any of phase 1,2,3,4 in the year up 
        to and including self._year
        We return factor times as many negative examples as we have positive examples
        """
        n_pos_examples = len(positive_df)
        n_neg_examples = n_pos_examples * factor
        kinase_list = [ geneid for _, geneid in self._symbol_to_id_map.items()]
        cancer_id_list = self._mesh_list
        positive_links = Link.fromDataFrameToLinkSet(positive_df)
        negative_links = set()
        i = 0  # use i to limit the number of attempts in case there is some problem
        while len(negative_links) < n_neg_examples and i < 1e6:
            i += 1
            random_cancer = random.choice(cancer_id_list)
            random_kinase = random.choice(kinase_list)
            randomLink = Link(kinase=random_kinase, cancer=random_cancer)
            if randomLink in positive_links:
                print("Skipping random link(%s,%s) since we found it in the positive set" %(random_kinase,random_cancer) )
                continue
            if randomLink in negative_links:
                print("Skipping random link (%s,%s) since we already added it to the negative set" %(random_kinase,random_cancer))
                continue
            negative_links.add(randomLink)
        print("[INFO] We generated a negative set with %d examples (the positive set has %d)" % (len(negative_links), len(positive_links)))
        # format as a pandas dataframe
        negative_dict_list = [ l.to_dict() for l in negative_links]
        return pd.DataFrame(negative_dict_list)

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
