import random

from .kcet_parser import KcetParser, PkPkiFilter
from .ct_by_phase_parser import CTParserByPhase

import pandas as pd
import numpy as np
import datetime
from typing import Set, Tuple, List
import os
import logging

logging.basicConfig(format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
                    datefmt='%Y-%m-%d:%H:%M:%S',
                    filename='kcet.log',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)


class Link:
    """
    Simple class that is intended for use to keep track of positive and negative links using a Hash
    For our use case, the _kinase field is NCBI Gene IDs for protein kinases and the _cancer field is
    a MeSH Id for a cancer (descendant of neoplasms)
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

    @staticmethod
    def fromEmbeddingsToLinkSet(df: pd.DataFrame) -> Set:
        """
        This method goes from the indices of the concept embeddings to Link objects.
        """
        linkset = set()
        for i, row in df.iterrows():
            # i (the index) is like this ncbigene7010-meshd018195
            k, c = i.split("-")
            L = Link(kinase=k, cancer=c)
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
    Functions such as get_training_and_test_embeddings refer to all phases; get_training_and_test_emebddings_phase4 is restricted to phase 4.
    """

    def __init__(self, clinical_trials: str, embeddings: str, words: str, n_pk: int = 5) -> None:
        kcetParser = KcetParser()
        self._pki_to_kinase_df = kcetParser._get_pki_to_kinase_list_dict_max_pk(n_pk=n_pk)
        if not isinstance(self._pki_to_kinase_df, pd.DataFrame):
            raise ValueError("_pki_to_kinase_dict needs to be a DataFrame")
        self._symbol_to_id_map = kcetParser.get_symbol_to_id_map()
        self._mesh_list = kcetParser.get_mesh_id_list()
        parser = CTParserByPhase(clinical_trials=clinical_trials)
        self._df_allphases = parser.get_all_phases(removeRedundantEntries=True)  # all positive data, phase 1,2,3,4
        self._df_phase4 = parser.get_phase_4(removeRedundantEntries=True)  # all positive data, phase 4 only
        self._n_pk = n_pk
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
        logger.info(
            "We ingested %d labeled word vectors from %s and %s" % (len(self._embeddings_df), embeddings, words))
        self._ncbigene2symbol_map = kcetParser.get_id_to_symbol_map()
        logger.info("We ingested %d symbol/NCBI gene id mappings" % (len(self._ncbigene2symbol_map)))
        self._meshid2disease_map = kcetParser.get_mesh_to_disease_map()
        logger.info("We ingested %d meshId/disease mapping" % (len(self._meshid2disease_map)))

    def get_words(self):
        return self._embeddings_df.index

    def get_embeddings(self) -> pd.DataFrame:
        """
        return a Pandas dataframe whose index is the words, and whose columns are the dimensions of the embeddings
        """
        return self._embeddings_df

    def get_training_and_test_embeddings(self, target_year: int, begin_year: int, end_year: int, factor: int = 10) -> \
            Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Get positive and negative training data for the target year
        Get data for test that goes from begin_year to end_year
        Return the data inform of pandas dataframes with word/concept embeddings representing the 'difference vectors'
        (see manuscript for details and definitions).
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

        positive_training_df = self.get_pos_training_embeddings(target_year=target_year)
        n_pos_train = len(positive_training_df)
        n_neg_train = n_pos_train * factor
        negative_training_df = self.get_neg_training_embeddings(target_year=target_year, n_neg_examples=n_neg_train)
        pos_test = self.get_positive_test_embeddings(target_year=target_year, begin_year=begin_year, end_year=end_year)
        n_neg_test = factor * len(pos_test)
        neg_test = self.get_negative_test_embeddings(negative_df=negative_training_df, year=target_year,
                                                     n_negative_test=n_neg_test)
        return positive_training_df, negative_training_df, pos_test, neg_test

    def get_training_and_test_embeddings_phase_4(self, target_year: int, begin_year: int, end_year: int,
                                                 factor: int = 10) -> \
            Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        The function is analogous to get_training_and_test_data (see this for documehtation) but is limited to phase 4 clinical studies
        """
        positive_training_df = self.get_pos_training_embeddings(target_year=target_year)
        n_pos_train = len(positive_training_df)
        n_neg_train = n_pos_train * factor
        negative_training_df = self.get_neg_training_embeddings(target_year=target_year, n_neg_examples=n_neg_train)
        if end_year < begin_year:
            raise ValueError("End year cannot be before start year")
        if begin_year < target_year:
            raise ValueError("Begin year cannot be before target year")
        positive_test_df = self.get_positive_test_embeddings(target_year=target_year, begin_year=begin_year,
                                                             end_year=end_year, phase4=True)
        n_pos_test = factor * len(positive_test_df)
        negative_test_df = self.get_negative_test_embeddings(negative_df=negative_training_df, year=target_year,
                                                             n_negative_test=n_pos_test)

        return positive_training_df, negative_training_df, positive_test_df, negative_test_df

    def get_pos_training_embeddings(self, target_year: int) -> pd.DataFrame:
        """
        get the positive embeddings for training
        """
        pos_train_df = self._get_positive_training_data_set(year=target_year)
        pos_train_vectors = self.get_disease_kinase_difference_vectors(pos_train_df)
        return pos_train_vectors

    def get_neg_training_embeddings(self, target_year: int, n_neg_examples: int) -> pd.DataFrame:
        """
        get negative embeddings for training
        We do so by choosing from among embeddings that are not positive
        n_neg_examples: number of embeddings to get (in general, 10 times the number of positive data for training)
        We take Random non-links that were not listed in any of phase 1,2,3,4 in the year up
        to and including self._year
        """
        kinase_list = [geneid for _, geneid in self._symbol_to_id_map.items()]
        cancer_id_list = self._mesh_list
        # The following help to keep track of positive examples
        positive_links = Link.fromDataFrameToLinkSet(self._df_allphases[self._df_allphases['year'] <= target_year])
        negative_links = set()
        n_skipped_link = 0
        i = 0  # use i to limit the number of attempts in case there is some problem
        df = pd.DataFrame(columns=self._embeddings_df.columns)
        while len(df) < n_neg_examples and i < 1e6:
            i += 1
            mesh_id = random.choice(cancer_id_list)
            ncbigene_id = random.choice(kinase_list)
            randomLink = Link(kinase=ncbigene_id, cancer=mesh_id)
            # Do not add a link from the positive set to the negative set 
            # Note that this can happen by chance and is not worrisome, but we log it
            if randomLink in positive_links or randomLink in negative_links:
                n_skipped_link += 1
                continue
            negative_links.add(randomLink)
            ncbigene_id_embedding = None
            mesh_id_embedding = None
            if ncbigene_id in self._embeddings_df.index:
                ncbigene_id_embedding = self._embeddings_df.loc[ncbigene_id]
            if mesh_id in self._embeddings_df.index:
                mesh_id_embedding = self._embeddings_df.loc[mesh_id]
            if ncbigene_id_embedding is not None and mesh_id_embedding is not None:
                diff_kinase_mesh = np.subtract(ncbigene_id_embedding, mesh_id_embedding)
                label = "%s-%s" % (ncbigene_id, mesh_id)
                df.loc[label] = diff_kinase_mesh
            i += 1
            if i % 10000 == 0 and i > 0:
                logger.info(
                    "Created %d/%d (%.1f%%) difference vectors" % (i, n_neg_examples, 100.0 * i / n_neg_examples))
        logger.info("Extracted %s kinase-cancer difference vectors" % len(df))
        return df

    def get_positive_test_embeddings(self, target_year: int, begin_year: int, end_year: int,
                                     phase4: bool = False) -> pd.DataFrame:
        """
        Get all of the positive examples from begin_year to end_year (inclusive).
        -- used for test in historical experiments
        """

        if phase4:
            within_valid_year_range = (self._df_phase4['year'] >= begin_year) & (self._df_phase4['year'] <= end_year)
            df_pos_test = self._df_phase4[within_valid_year_range]
        else:
            within_valid_year_range = (self._df_allphases['year'] >= begin_year) & (
                    self._df_allphases['year'] <= end_year)
            df_pos_test = self._df_allphases[within_valid_year_range]
        all_positive_links_up_to_target = self.get_all_phases_all_pk_pki(target_year=target_year)
        positive_test_links = set()
        n_skipped = 0
        for _, row in df_pos_test.iterrows():
            cancer_mesh_id = row['mesh_id']
            kinase_ncbi_gene_id = row['gene_id']
            candidate = Link(cancer=cancer_mesh_id, kinase=kinase_ncbi_gene_id)
            if candidate not in all_positive_links_up_to_target:
                positive_test_links.add(candidate)
            else:
                n_skipped += 1
        logger.info("{} candidate PK/cancer pairs were skipped for the positive test set".format(n_skipped))
        kinase_list = []
        cancer_list = []
        n_skipped_link = 0
        df = pd.DataFrame(columns=self._embeddings_df.columns)
        for link in positive_test_links:
            if link in all_positive_links_up_to_target:
                # do not include a positive example if it was already known at training time!
                n_skipped_link += 1
                continue
            ncbigene_id = link.kinase
            mesh_id = link.cancer
            kinase_list.append(ncbigene_id)
            cancer_list.append(mesh_id)
            ncbigene_id_embedding = None
            mesh_id_embedding = None
            if ncbigene_id in self._embeddings_df.index:
                ncbigene_id_embedding = self._embeddings_df.loc[ncbigene_id]
            if mesh_id in self._embeddings_df.index:
                mesh_id_embedding = self._embeddings_df.loc[mesh_id]
            if ncbigene_id_embedding is not None and mesh_id_embedding is not None:
                diff_kinase_mesh = np.subtract(ncbigene_id_embedding, mesh_id_embedding)
                label = "%s-%s" % (ncbigene_id, mesh_id)
                df.loc[label] = diff_kinase_mesh
        logger.info(
            "Skipped %d links for testing that were already present in training data (expected behavior)" % n_skipped_link)
        return df

    def get_negative_test_embeddings(self, negative_df: pd.DataFrame, year: int, n_negative_test) -> pd.DataFrame:
        """
        Get negative examples after the target year-- used for testing
        in historical experiments. Note that we take examples that are negative
        from the perspective of the current time -- we are taking factor-times more negative
        examples than positive examples, and this function chooses a set that is distinct
        from the set of examples use prior to the target year (negative_df).
        """

        kinase_list = [geneid for _, geneid in self._symbol_to_id_map.items()]
        cancer_id_list = self._mesh_list
        positive_links = Link.fromDataFrameToLinkSet(self._df_allphases[self._df_allphases['year'] <= year])
        pretarget_negative_links = Link.fromEmbeddingsToLinkSet(negative_df)
        negative_links = set()
        n_skipped_link = 0
        df = pd.DataFrame(columns=self._embeddings_df.columns)
        i = 0  # use i to limit the number of attempts in case there is some problem
        while len(df) < n_negative_test and i < 1e6:
            i += 1
            mesh_id = random.choice(cancer_id_list)
            ncbigene_id = random.choice(kinase_list)
            randomLink = Link(kinase=ncbigene_id, cancer=mesh_id)
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
            ncbigene_id_embedding = None
            mesh_id_embedding = None
            if ncbigene_id in self._embeddings_df.index:
                ncbigene_id_embedding = self._embeddings_df.loc[ncbigene_id]
            if mesh_id in self._embeddings_df.index:
                mesh_id_embedding = self._embeddings_df.loc[mesh_id]
            if ncbigene_id_embedding is not None and mesh_id_embedding is not None:
                diff_kinase_mesh = np.subtract(ncbigene_id_embedding, mesh_id_embedding)
                label = "%s-%s" % (ncbigene_id, mesh_id)
                df.loc[label] = diff_kinase_mesh

        logger.info("Skipped %d links that were found previously (expected behavior)" % n_skipped_link)
        logger.info("We generated a negative test set with %d examples (the positive set has %d)" % (
            len(negative_links), len(positive_links)))
        return df

    def get_all_phases_all_pk_pki(self, target_year: int):
        """
        It is a conservative estimate to assume that all connections between a PK and PK are valid for the testing set
        In contrast, for training, we use only validated items (phase 4)
        This method returns a data frame with cancer/PK links that are derived from all studies up to the target year,
        and all PK/PKI links (i.e., not limited to the n_pk_pki parameter that is used for the training set)
        """
        pkpki = PkPkiFilter()
        all_pk_pki_df = pkpki.get_all_pk_pki()
        # all_pk_pki_df has rows like this - PKI:abemaciclib; PK: CDK4, ACT_VALUE:0.0000599..., PMID:24919854
        all_links = set()
        all_phases = self._df_allphases[self._df_allphases['year'] <= target_year]
        for _, row in all_phases.iterrows():
            kinase = row['kinase']  # e.g., CDK4
            gene_id = row['gene_id']
            cancer = row['mesh_id']
            pk_pki = all_pk_pki_df[all_pk_pki_df['PK'] == kinase]
            for idx, item in pk_pki.iterrows():
                pk = item['PK']
                link = Link(cancer=cancer, kinase=gene_id)
                all_links.add(link)
        return all_links

    def get_data_for_novel_prediction(self, target_year: int, factor: int = 10) -> Tuple[
        pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        This method creates positive and negative training sets including everything up to the current year
        It then creates examples for all other predictions that we will use for the novel predictions
        It returns three dataframes with embeddings.
        """
        positive_training_df = self.get_pos_training_embeddings(target_year=target_year)
        n_neg = factor * len(positive_training_df)
        negative_training_df = self.get_neg_training_embeddings(target_year=target_year, n_neg_examples=n_neg)
        prediction_df = pd.DataFrame(columns=self._embeddings_df.columns)
        kinase_list = [geneid for _, geneid in self._symbol_to_id_map.items()]
        cancer_id_list = self._mesh_list
        i = 0
        total = len(kinase_list) * len(cancer_id_list)
        print("Links to be extracted: {}".format(total))
        # We remove all positive protein-kinase/cancer associations regardless of phase
        positive_links = self.get_all_phases_all_pk_pki(target_year=target_year)
        negative_links = Link.fromEmbeddingsToLinkSet(negative_training_df)
        for ncbigene_id in kinase_list:
            for mesh_id in cancer_id_list:
                i += 1
                if i % 2000 == 0:
                    print("\r{}/{} links extracted ({:.2f}%)".format(i, total, 100 * i / total), end="")
                randomLink = Link(kinase=ncbigene_id, cancer=mesh_id)
                if randomLink in positive_links or randomLink in negative_links:
                    continue
                ncbigene_id_embedding = None
                mesh_id_embedding = None
                if ncbigene_id in self._embeddings_df.index:
                    ncbigene_id_embedding = self._embeddings_df.loc[ncbigene_id]
                if mesh_id in self._embeddings_df.index:
                    mesh_id_embedding = self._embeddings_df.loc[mesh_id]
                if ncbigene_id_embedding is not None and mesh_id_embedding is not None:
                    diff_kinase_mesh = np.subtract(ncbigene_id_embedding, mesh_id_embedding)
                    label = "%s-%s" % (ncbigene_id, mesh_id)
                    prediction_df.loc[label] = diff_kinase_mesh
        return positive_training_df, negative_training_df, prediction_df

    def _get_positive_training_data_set(self, year: int) -> pd.DataFrame:
        """
        Positive training set: all links of phase 4 up to the year given in the constructor
        """
        if not isinstance(year, int):
            raise ValueError("year must be an integer")
        return self._df_phase4[self._df_phase4['year'] <= year]

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
        if "gene_id" not in examples.columns:
            raise ValueError("Input dataframe must contain a column called gene_id")
        if "mesh_id" not in examples.columns:
            raise ValueError("Input dataframe must contain a column called mesh_id")
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
                label = "%s-%s" % (ncbigene_id, mesh_id)
                df.loc[label] = diff_kinase_mesh
            i += 1
            if i % 10000 == 0 and i > 0:
                logger.info("Created %d/%d (%.1f%%) difference vectors" % (i, total, 100.0 * i / total))
        logger.info("Extracted %s kinase-cancer difference vectors" % len(df))
        logger.info("Initial data: %d examples" % len(examples))
        logger.info("Could not identify %d gene ids" % len(unidentified_genes))
        logger.info("Could not identify %d MeSH ids" % len(unidentified_cancers))
        return df

    def get_summary(self) -> pd.DataFrame:
        """
        Return a data frame with counts and descriptive statistics about the current dataset that can be used
        for display in a Jupyter notebook etc.
        """
        data = []
        n_pki = len(self._pki_to_kinase_df)
        data.append(['Above threshold PKI/PKI links', "{:d}".format(n_pki)])
        df = self._pki_to_kinase_df.groupby("PKI")["PK"].count()
        mean_pk_per_pki = np.mean(df)
        data.append(['mean PKs per PKI (DrugCentral dataset)', "{:.2f}".format(mean_pk_per_pki)])
        n_unique_kinases = len(pd.unique(self._pki_to_kinase_df['PK']))
        data.append(['protein kinases in DrugCentral dataset', "{:d}".format(n_unique_kinases)])
        n_unique_PKIs = len(pd.unique(self._pki_to_kinase_df['PKI']))
        data.append(['protein kinase inhibitors in DrugCentral dataset', "{:d}".format(n_unique_PKIs)])
        n_ncbi_kinases = len(self._symbol_to_id_map)
        data.append(['protein kinases in NCBI gene dataset', "{:d}".format(n_ncbi_kinases)])
        n_mesh = len(self._mesh_list)
        data.append(['MeSH Ids for cancer concepts', "{:d}".format(n_mesh)])
        n_all_phases = len(self._df_allphases)
        data.append(['Clinical trials included in this study', "{:d}".format(n_all_phases)])
        for i in [1, 2, 3, 4]:
            phase = 'Phase {}'.format(i)
            n = len(self._df_allphases[self._df_allphases['phase'] == phase])
            message = 'Phase {} trials included in this analysis'.format(i)
            data.append([message, "{:d}".format(n)])
        n_embeddings = len(self._embeddings_df)
        data.append(['word/concept embeddings', "{:d}".format(n_embeddings)])
        return pd.DataFrame(data, columns=['Item', 'Value'])
