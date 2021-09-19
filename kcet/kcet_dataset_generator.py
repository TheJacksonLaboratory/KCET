import random
from .kcet_parser import KcetParser
from .ct_by_phase_parser import CTParserByPhase

import pandas as pd
import numpy as np
import datetime
from typing import Set, Tuple
import os
import logging

logging.basicConfig(filename='kcet.log', level=logging.INFO)


class Link:
    """
    Simple class that is intended for use to keep track of positive and negative links using a Hash
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

    def __str__(self):
        return self._cancer + "-" + self._kinase

    @staticmethod
    def fromDataFrameToLinkSet(df: pd.DataFrame) -> Set:
        linkset = set()
        for index, row in df.iterrows():
            m = row['mesh_id']
            g = row['gene_id']
            L = Link(kinase=g, cancer=m)
            linkset.add(L)
        return linkset


def get_current_year():
    now = datetime.datetime.now()  # default to current year
    year = now.year
    return year


class KcetDatasetGenerator:
    """
    Class to generate test, training, and prediction files.
    Use of terminology
    – Training set: A set of examples used for learning.
    – Validation set: A set of examples used to tune the parameters of a classifier.
    – Test set: A set of examples used only to assess the performance of a fully-specified classifier.
    In our case, we use historical data for training/test by extracting word/concept embedddings up to the target year (e.g., 2010)
    and use subsequent years for testing (e.g., data from after 2010).
    For some experiments, we restrict our analysis to phase 4 clinical trials. For others, we take all phases.
    Functions such as get_training_and_test_data refer to all phases; get_training_and_test_data_phase4 is restricted to phase 4.
    """

    def __init__(self, clinical_trials: str, embeddings: str, words: str) -> None:
        kcetParser = KcetParser()
        self._pki_to_kinase_dict = kcetParser.get_pki_to_kinase_list_dict()
        self._symbol_to_id_map = kcetParser.get_symbol_to_id_map()
        self._mesh_list = kcetParser.get_mesh_id_list()
        parser = CTParserByPhase(clinical_trials=clinical_trials)
        self._df_allphases = parser.get_all_phases(removeRedundantEntries=True)  # all positive data, phase 1,2,3,4
        self._df_phase4 = parser.get_phase_4(removeRedundantEntries=True)  # all positive data, phase 4 only
        # add the embeddings
        if not os.path.exists(embeddings):
            raise FileNotFoundError("Could not find embedding file at %s" % embeddings)
        if not os.path.exists(words):
            raise FileNotFoundError("Could not find words file at %s" % words)
        embedding = np.load(embeddings, mmap_mode=None, allow_pickle=False, fix_imports=True, encoding='ASCII')
        word_list = []
        with open(words) as f:
            for line in f:
                word = line[2:-3]
                word_list.append(word)
        self._embeddings_df = pd.DataFrame(data=embedding, index=word_list)
        logging.info(
            "We ingested %d labeled word vectors from %s and %s" % (len(self._embeddings_df), embeddings, words))
        self._ncbigene2symbol_map = kcetParser.get_id_to_symbol_map()
        logging.info("We ingested %d symbol/NCBI gene id mappings" % (len(self._ncbigene2symbol_map)))
        self._meshid2disease_map = kcetParser.get_mesh_to_disease_map()
        logging.info("We ingested %d meshId/disease mapping" % (len(self._meshid2disease_map)))

    def get_words(self):
        return self._embeddings_df.index

    def get_embeddings(self) -> pd.DataFrame:
        """
        return a Pandas dataframe whose index is the words, and whose columns are the dimensions of the embeddings
        """
        return self._embeddings_df

    def get_training_and_test_data(self, target_year: int, begin_year: int, end_year: int, factor: int = 10) -> \
            Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        WAS get_data_years_after_target_year_upto_later_year
        Get positive and negative training data for the target year 
        Get data for test that goes from begin_year to end_year
        Parameters
        ----------
        target_year : int
            The year up to which embeddings were created from PubMed abstracts. Training is done with clinical trials data up to the target year
        begin_year : int
            The first year of the test data
        end_year : int
            The last year of the test data
        factor : int
            We randomly choose factor-times as many negative training examples as there is positive training examples
        """
        positive_training_df = self._get_positive_training_data_set(year=target_year)
        negative_training_df = self._get_negative_training_dataset(pos_training_df=positive_training_df,
                                                                   year=target_year, factor=factor)
        if end_year < begin_year:
            raise ValueError("End year cannot be before start year")
        if begin_year < target_year:
            raise ValueError("Begin year cannot be before target year")

        positive_test_df = self._get_positive_test_data(target_year=target_year, begin_year=begin_year,
                                                        end_year=end_year)
        negative_test_df = self._get_negative_test_data(negative_df=negative_training_df, year=target_year)
        return positive_training_df, negative_training_df, positive_test_df, negative_test_df

    def get_data_years_after_target_year_upto_later_year_phase_4(self, target_year: int, mid_year: int,
                                                                 factor: int = 10,
                                                                 num_years_later: int = 1) -> Tuple[pd.DataFrame]:
        raise Exception("Not supported")

    def get_training_and_test_data_phase_4(self, target_year: int, begin_year: int, end_year: int, factor: int = 10) -> \
            Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        The function is analogous to get_training_and_test_data (see this for documehtation) but is limited to phase 4 clinical studies
        """
        pos_train_df = self._get_positive_training_data_set(year=target_year)
        neg_train_df = self._get_negative_training_dataset(pos_training_df=pos_train_df, year=target_year,
                                                           factor=factor)
        if end_year < begin_year:
            raise ValueError("End year cannot be before start year")
        if begin_year < target_year:
            raise ValueError("Begin year cannot be before target year")
        positive_test_df = self._get_positive_test_data_phase_4(target_year=target_year, begin_year=begin_year,
                                                                end_year=end_year)
        negative_test_df = self._get_negative_test_data(negative_df=neg_train_df, year=target_year)
        return pos_train_df, neg_train_df, positive_test_df, negative_test_df

    def get_data_for_novel_prediction(self, current_year: int, factor: int = 10):
        if current_year != get_current_year():
            raise Exception("For novel prediction, year must be the current year! ")
        positive_training_df = self._get_positive_training_data_set(year=current_year)
        negative_training_df = self._get_negative_training_dataset(pos_training_df=positive_training_df,
                                                                   year=current_year, factor=factor)
        prediction_df = self.get_prediction_dataset()
        return positive_training_df, negative_training_df, prediction_df

    def _get_positive_training_data_set(self, year: int) -> pd.DataFrame:
        """
        Positive training set: all links of phase 4 up to the year given in the constructor
        """
        if not isinstance(year, int):
            raise ValueError("year must be an integer")
        return self._df_phase4[self._df_phase4['year'] <= year]

    def _get_positive_test_data(self, target_year: int, begin_year: int, end_year: int) -> pd.DataFrame:
        """
        Get all of the positive examples from begin_year to end_year (inclusive).
        -- used for test in historical experiments
        """
        within_valid_year_range = (self._df_allphases['year'] >= begin_year) & (self._df_allphases['year'] <= end_year)
        df_pos_test = self._df_allphases[within_valid_year_range]
        positive_test_links = Link.fromDataFrameToLinkSet(df_pos_test)
        # Ground-truth training data is up to the target year only!
        all_phases_positive_links = Link.fromDataFrameToLinkSet(
            self._df_allphases[self._df_allphases['year'] <= target_year])

        kinase_list = []
        cancer_list = []
        n_skipped_link = 0
        for link in positive_test_links:
            if link in all_phases_positive_links:
                # do not include a positive example if it was already known at training time!
                n_skipped_link += 1
            else:
                kinase = link.kinase
                cancer = link.cancer
                kinase_list.append(kinase)
                cancer_list.append(cancer)
        logging.info(
            "Skipped %d links for testing that were already present in training data (expected behavior)" % n_skipped_link)
        return pd.DataFrame(list(zip(cancer_list, kinase_list)), columns=['mesh_id', 'gene_id'])

    def _get_positive_test_data_phase_4(self, target_year: int, begin_year: int, end_year: int) -> pd.DataFrame:
        """
        Get all of the positive links (of phase 4)  after the target year-- used for test
        in historical experiments
        """
        within_valid_year_range = (self._df_phase4['year'] >= begin_year) & (self._df_phase4['year'] <= end_year)
        df_pos_test = self._df_phase4[within_valid_year_range]
        positive_test_links = Link.fromDataFrameToLinkSet(df_pos_test)
        all_phases_positive_links = Link.fromDataFrameToLinkSet(
            self._df_allphases[self._df_allphases['year'] <= target_year])

        kinase_list = []
        cancer_list = []
        n_skipped_link = 0
        for link in positive_test_links:
            if link in all_phases_positive_links:
                n_skipped_link += 1
            else:
                kinase = link.kinase
                cancer = link.cancer
                kinase_list.append(kinase)
                cancer_list.append(cancer)
        logging.info(
            "Skipped %d links for testing that were already present in training data (expected behavior)" % n_skipped_link)
        return pd.DataFrame(list(zip(cancer_list, kinase_list)), columns=['mesh_id', 'gene_id'])

    def _get_negative_test_data(self, negative_df: pd.DataFrame, year: int) -> pd.DataFrame:
        """
        Get negative examples after the target year-- used for testing
        in historical experiments. Note that we take examples that are negative
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
        n_skipped_link = 0
        i = 0  # use i to limit the number of attempts in case there is some problem
        while len(negative_links) < n_neg_examples and i < 1e6:
            i += 1
            random_cancer = random.choice(cancer_id_list)
            random_kinase = random.choice(kinase_list)
            randomLink = Link(kinase=random_kinase, cancer=random_cancer)
            if randomLink in positive_links:
                n_skipped_link += 1
                continue
            if randomLink in pretarget_negative_links:
                n_skipped_link += 1
                continue
            if randomLink in negative_links:
                n_skipped_link += 1
                continue
            negative_links.add(randomLink)
        logging.info("Skipped %d links that were found previously (expected behavior)" % n_skipped_link)
        logging.info("We generated a negative test set with %d examples (the positive set has %d)" % (
            len(negative_links), len(positive_links)))
        # format as a pandas dataframe
        negative_dict_list = [l.to_dict() for l in negative_links]
        return pd.DataFrame(negative_dict_list)

    def _get_unique_positive_training_count(self, pos_training_df: pd.DataFrame) -> int:
        """
        The pandas dataframe with positive training examples can have duplicates if we are 
        training for phases 1,2,3,4 (all phases are combined into cancer-protein kinase pairs in further processing).
        In order to get the right number of negative training examples, we need to find the number
        of unique positive examples, which is the purpose of this function.
        We do this simply by making a key out of the mesh and gene ids and keeping track of this.
        """
        if not isinstance(pos_training_df, pd.DataFrame):
            raise ValueError("pos training df must be a pandas dataframe")
        unique_pairs = set()
        for i, row in pos_training_df.iterrows():
            key = row['mesh_id'] + '-' + row['gene_id']
            unique_pairs.add(key)
        return len(unique_pairs)

    def _get_negative_training_dataset(self, pos_training_df: pd.DataFrame, year: int,
                                       factor: int = 10) -> pd.DataFrame:
        """
        Get a negative training set.
        We take Random non-links that were not listed in any of phase 1,2,3,4 in the year up 
        to and including self._year
        We return factor times as many negative examples as we have positive examples. Positive examples are links from phase 4 only!
        """

        kinase_list = [geneid for _, geneid in self._symbol_to_id_map.items()]
        cancer_id_list = self._mesh_list
        positive_links = Link.fromDataFrameToLinkSet(self._df_allphases[self._df_allphases['year'] <= year])
        # positive_training_links = Link.fromDataFrameToLinkSet(pos_training_df)
        n_pos_examples = self._get_unique_positive_training_count(
            pos_training_df)  # number of links in positive training set
        n_neg_examples = n_pos_examples * factor
        negative_links = set()
        n_skipped_link = 0
        i = 0  # use i to limit the number of attempts in case there is some problem
        while len(negative_links) < n_neg_examples and i < 1e6:
            i += 1
            random_cancer = random.choice(cancer_id_list)
            random_kinase = random.choice(kinase_list)
            randomLink = Link(kinase=random_kinase, cancer=random_cancer)
            # Do not add a link from the positive set to the negative set 
            # Note that this can happen by chance and is not worrisome, but we log it
            if randomLink in positive_links or randomLink in negative_links:
                n_skipped_link += 1
                continue
            negative_links.add(randomLink)
        logging.info(
            "Skipped %d random links that were found in previously defined sets (expected behavior)" % n_skipped_link)
        logging.info("[INFO] We generated a negative training set with %d examples (the positive set has %d)" % (
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
        df = pd.DataFrame(columns=self._embeddings_df.columns)
        total = len(examples.index)
        if total == 0:
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
                # diff_kinase_mesh_list_pos_train.append(diff_kinase_mesh)
                # diff_index_pos_train.append(ncbigene_id + "," + mesh_id)
                label = "%s-%s" % (ncbigene_id, mesh_id)
                df.loc[label] = diff_kinase_mesh
            i += 1
            if i % 10000 == 0 and i > 0:
                logging.info("Created %d/%d (%.1f%%) difference vectors" % (i, total, 100.0 * i / total))
        logging.info("Extracted %s kinase-cancer difference vectors" % len(df))
        logging.info("Initial data: %d examples" % len(examples))
        logging.info("Could not identify %d gene ids" % len(unidentified_genes))
        logging.info("Could not identify %d MeSH ids" % len(unidentified_cancers))
        return df
