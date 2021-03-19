import random
import pandas as pd
import datetime
from typing import List, Tuple
import logging
logging.basicConfig(filename='kcet.log', level=logging.INFO)

from .kcet_parser import KcetParser
from .ct_by_phase_parser import CTParserByPhase


class Link:
    """
    Simple class to keep track of positive and negative links using a Hash
    """

    def __init__(self, kinase: str, cancer: str) -> None:
        self._kinase = kinase
        self._cancer = cancer

    def __hash__(self):
        return hash((self._kinase, self._cancer))

    def __eq__(self, other):
        return self.__class__ == other.__class__ and self._kinase == other.kinase and self._cancer == other.cancer

    @property
    def kinase(self):
        return self._kinase

    @property
    def cancer(self):
        return self._cancer

    def to_dict(self):
        return {'mesh_id': self._cancer, 'gene_id': self._kinase}

    @staticmethod
    def fromDataFrameToLinkSet(df: pd.DataFrame) -> List:
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
        kcetParser = KcetParser()
        self._pki_to_kinase_dict = kcetParser.get_pki_to_kinase_list_dict()
        self._symbol_to_id_map = kcetParser.get_symbol_to_id_map()
        self._mesh_list = kcetParser.get_mesh_id_list()
        parser = CTParserByPhase(clinical_trials=clinical_trials)
        self._df_allphases = parser.get_all_phases()  # all positive data, phase 1,2,3,4
        self._df_phase4 = parser.get_phase_4()  # all positive data, phase 4 only

    def get_current_year(self):
        now = datetime.datetime.now()  # default to current year
        year = now.year
        return year

    def get_data_for_target_year(self, target_year: int, factor: int = 10) -> Tuple[pd.DataFrame, pd.DataFrame,
                                                                                    pd.DataFrame, pd.DataFrame]:
        """
        Get positive and negative data for the target year indicated by the argument.
        """
        positive_df = self._get_positive_data_set(year=target_year)
        negative_df = self._get_negative_training_dataset(pos_training_df = positive_df, year=target_year, factor=factor)
        current_year = self.get_current_year()
        if target_year < current_year:  # historical prediction
            positive_validation_df = self._get_positive_validation_data_set(year=target_year)
            negative_validation_df = self._get_negative_validation_data_set(negative_df=negative_df, year=target_year)
            return positive_df, negative_df, positive_validation_df, negative_validation_df
        else:  # novel prediction
            prediction_df = self.get_prediction_dataset()
            return positive_df, negative_df, prediction_df

    def get_data_for_target_year_and_later_year(self, target_year: int, factor: int = 10, num_years_later: int = 1) -> \
            Tuple[pd.DataFrame, pd.DataFrame,pd.DataFrame, pd.DataFrame]:
        """
        Get positive and negative training for the target year and positive and negative validation for num_years_later years later than the target year.
        """
        positive_df = self._get_positive_data_set(year=target_year)
        negative_df = self._get_negative_training_dataset(pos_training_df = positive_df,year=target_year, factor=factor)
        current_year = self.get_current_year()
        new_target_year = target_year + num_years_later
        if target_year < current_year and new_target_year < current_year:  # historical prediction
            positive_validation_df = self._get_positive_validation_data_set_later_year(old_target_year=target_year,new_target_year=new_target_year)
            negative_validation_df = self._get_negative_validation_data_set(negative_df=negative_df,year=target_year)
            return positive_df, negative_df, positive_validation_df, negative_validation_df
        else:
            raise Exception("target year and the future year must be previous to the current year! ")

    def get_data_for_target_year_phase_4(self, target_year: int, factor: int = 10) -> Tuple[pd.DataFrame]:
        """
        Get positive and negative data for the target year indicated by the argument. The positive training and positive validation are
        links from phase 4.
        """
        positive_df = self._get_positive_data_set(year=target_year)
        negative_df = self._get_negative_training_dataset(pos_training_df = positive_df,year=target_year, factor=factor)
        now = datetime.datetime.now()  # default to current year
        currentyear = now.year
        if target_year < currentyear:  # historical prediction
            positive_validation_df = self._get_positive_validation_data_set_phase_4(year=target_year)  # get phase 4 for the years after the target year
            negative_validation_df = self._get_negative_validation_data_set(negative_df=negative_df, year=target_year)
            return positive_df, negative_df, positive_validation_df, negative_validation_df
        else:  # novel prediction
            prediction_df = self.get_prediction_dataset()
            return positive_df, negative_df, prediction_df

    def _get_positive_data_set(self, year: int) -> pd.DataFrame:
        """
        Positive training set: all links of phase 4 up to the year given in the constructor
        """
        if not isinstance(year, int):
            raise ValueError("year must be an integer")

        not_later_than_target_year = self._df_phase4['year'] <= year
        return self._df_phase4[not_later_than_target_year]

    def _get_positive_validation_data_set_phase_4(self, year: int) -> pd.DataFrame:
        """
        Get all of the positive links (of phase 4)  after the target year-- used for validation
        in historical validation experiments
        """
        print("GPVDS year=" + str(year))
        phase_4_later_than_target_year = self._df_phase4['year'] > year
        return self._df_phase4[phase_4_later_than_target_year]

    def _get_positive_validation_data_set(self, year: int) -> pd.DataFrame:
        """
        Get all of the positive examples from all phases after the target year-- used for validation
        in historical validation experiments
        """
        print("GPVDS year=" + str(year))
        later_than_target_year = self._df_allphases['year'] > year
        df_pos_validation = self._df_allphases[later_than_target_year]
        positive_validation_links = Link.fromDataFrameToLinkSet(df_pos_validation)
        kinase_list = []
        cancer_list = []
        for link in positive_validation_links:
             kinase = link.kinase
             cancer = link.cancer
             kinase_list.append(kinase)
             cancer_list.append(cancer)
        df_pos_validation = pd.DataFrame(list(zip( cancer_list,kinase_list)), columns =['mesh_id', 'gene_id'])
        return df_pos_validation

    def _get_positive_validation_data_set_later_year(self, old_target_year: int, new_target_year: int) -> pd.DataFrame:
        """
        Get all of the positive examples from one year after the old_target_year until the new_target_year.
        For example, if the old_target_year is 2010 and new_target_year is 2015, then the positive validation set is all positive links from
        2011 until 2015.
        -- used for validation in historical validation experiments
        """
        data_until_new_target_year=(self._df_allphases['year']>old_target_year) & (self._df_allphases['year'] <= new_target_year)
        df_pos_validation = self._df_allphases[data_until_new_target_year]
        positive_validation_links = Link.fromDataFrameToLinkSet(df_pos_validation)
        kinase_list = []
        cancer_list = []
        for link in positive_validation_links:
            kinase = link.kinase
            cancer = link.cancer
            kinase_list.append(kinase)
            cancer_list.append(cancer)
        df_pos_validation = pd.DataFrame(list(zip(kinase_list, cancer_list)), columns=['kinase', 'cancer'])
        return df_pos_validation

    def _get_negative_validation_data_set(self,  negative_df: pd.DataFrame, year: int) -> pd.DataFrame:
        """
        Get negative examples after the target year-- used for validation
        in historical validation experiments. Note that we take examples that are negative
        from the perspective of the current time -- we are taking factor-times more negative
        examples than positive examples, and this function chooses a set that is distinct
        from the set of examples use prior to the target year (the same size as the
        negative_df that is passed to us).
        """

        kinase_list = [geneid for _, geneid in self._symbol_to_id_map.items()]
        cancer_id_list = self._mesh_list
        positive_links = Link.fromDataFrameToLinkSet(self._df_allphases[self._df_allphases['year'] <= year])
        n_neg_examples = len(negative_df)
        pretarget_negative_links = Link.fromDataFrameToLinkSet(negative_df)
        negative_links = set()
        i = 0  # use i to limit the number of attempts in case there is some problem
        while len(negative_links) < n_neg_examples and i < 1e6:
            i += 1
            random_cancer = random.choice(cancer_id_list)
            random_kinase = random.choice(kinase_list)
            randomLink = Link(kinase=random_kinase, cancer=random_cancer)
            if randomLink in positive_links:
                print("Skipping random link (%s,%s) since we found it in the positive set" % (
                    random_kinase, random_cancer))
                continue
            if randomLink in pretarget_negative_links:
                print("Skipping random link (%s,%s) since we already added it to the negative set" % (
                    random_kinase, random_cancer))
                continue
            negative_links.add(randomLink)
        print("[INFO] We generated a negative set with %d examples (the positive set has %d)" % (
            len(negative_links), len(positive_links)))
        # format as a pandas dataframe
        negative_dict_list = [l.to_dict() for l in negative_links]
        return pd.DataFrame(negative_dict_list)

    def _get_negative_training_dataset(self, pos_training_df:pd.DataFrame, year: int, factor: int = 10) -> pd.DataFrame:
        """
        Get a negative training set.
        We take Random non-links that were not listed in any of phase 1,2,3,4 in the year up 
        to and including self._year
        We return factor times as many negative examples as we have positive examples. Positive exsmples are links from phase 4 only!
        """

        kinase_list = [geneid for _, geneid in self._symbol_to_id_map.items()]
        cancer_id_list = self._mesh_list
        positive_links = Link.fromDataFrameToLinkSet(self._df_allphases[self._df_allphases['year'] <= year])
        positive_training_links = Link.fromDataFrameToLinkSet(pos_training_df)
        n_pos_examples = len(positive_training_links) # number of links in positive training set
        n_neg_examples = n_pos_examples * factor
        negative_links = set()
        i = 0  # use i to limit the number of attempts in case there is some problem
        while len(negative_links) < n_neg_examples and i < 1e6:
            i += 1
            random_cancer = random.choice(cancer_id_list)
            random_kinase = random.choice(kinase_list)
            randomLink = Link(kinase=random_kinase, cancer=random_cancer)
            if randomLink in positive_links: # check if the random link is in any of pahses 1,2,3,4 or not. If yes, generate another link
                print("Skipping random link(%s,%s) since we found it in the positive set" % (
                    random_kinase, random_cancer))
                continue
            if randomLink in negative_links:
                print("Skipping random link (%s,%s) since we already added it to the negative set" % (
                    random_kinase, random_cancer))
                continue
            negative_links.add(randomLink)
        print("[INFO] We generated a negative set with %d examples (the positive set has %d)" % (
            len(negative_links), len(positive_links)))
        # format as a pandas dataframe
        negative_dict_list = [l.to_dict() for l in negative_links]
        return pd.DataFrame(negative_dict_list)

    def get_prediction_dataset(self) -> pd.DataFrame:
        """
        This set contains all possible links between untargeted protein kinases
        and the cancers except the positive dataset calculated up to self._year
        """
        positive_links = Link.fromDataFrameToLinkSet(self._df_allphases)
        kinase_list = [geneid for _, geneid in self._symbol_to_id_map.items()]
        cancer_id_list = self._mesh_list
        prediction_list = []
        for kinase in kinase_list:
            for cancer in cancer_id_list:
                L = Link(cancer=cancer, kinase=kinase)
                if L in positive_links:
                    continue  # Do not include positive links in the prediction set
                prediction_list.append(L.to_dict())
        return pd.DataFrame(prediction_list)
