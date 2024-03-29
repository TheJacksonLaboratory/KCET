from .kcet_dataset_generator import KcetDatasetGenerator

import pandas as pd
import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV

import logging

logging.basicConfig(filename='kcet.log', level=logging.INFO)


class KcetRandomForest:
    """
    This class is a wrapper around scikit learn functions for random forest classification.
    """

    def __init__(self,
                 data_gen: KcetDatasetGenerator,
                 embedddingfile: str,
                 wordsfile: str,
                 target: int,
                 factor: int = 10) -> None:
        self._data_generator = data_gen
        self._target_year = target
        self._factor = factor
        if not os.path.isfile(embedddingfile):
            raise FileNotFoundError("Could not find embedding file at " + embedddingfile)
        if not os.path.isfile(wordsfile):
            raise FileNotFoundError("Could not find embedding/words file at " + wordsfile)

    def classify(self, begin_year: int, end_year: int, phase4: bool = False):
        """
        Perform random forest learning. From the vectors extracted from the data from the target year, predict
        clinical trials starting at midyear and going num_years_later
        For instance, 
        target_year = 2010
        mid_year = 2019
        num_years_after_the_mid_year = 1
        creates test tests from 2019 to 2020.
        """
        # The following four data frames contain the names of the cancers and protein kinases and other columns.
        #
        """
        if phase4:
            logging.info("classify-getting data for phase 4")
            pos_train_df, neg_train_df, pos_test_df, neg_test_df = \
                self._data_generator.get_training_and_test_data_phase_4(target_year=self._target_year, begin_year=begin_year,
                                                                end_year=end_year)
        else:
            logging.info("classify-getting data for all phases")
            pos_train_df, neg_train_df, pos_test_df, neg_test_df = \
                self._data_generator.get_training_and_test_data(target_year=self._target_year,
                                                                        begin_year=begin_year, end_year=end_year)
       
        # Now we need to extract the corresponding embedded vectors
        logging.info(
            "classify examples: pos train: {}, neg train {}, pos test {}, neg test {}"
                .format(len(pos_train_df), len(neg_train_df), len(pos_test_df), len(neg_test_df)))
        """
        if phase4:
            pos_train_vectors, neg_train_vectors, pos_test_vectors, neg_test_vectors = self._data_generator.get_training_and_test_embeddings_phase_4(
                self._target_year, begin_year=begin_year,
                end_year=end_year)
        else:
            pos_train_vectors, neg_train_vectors, pos_test_vectors, neg_test_vectors = self._data_generator.get_training_and_test_embeddings(
                self._target_year, begin_year=begin_year,
                end_year=end_year)

        n_pos_train = pos_train_vectors.shape[0]
        n_neg_train = neg_train_vectors.shape[0]
        X_train = pd.concat([pos_train_vectors, neg_train_vectors])
        y_train = np.concatenate((np.ones(n_pos_train), np.zeros(n_neg_train)))
        # Prepare for random forest testing
        n_pos_test = pos_test_vectors.shape[0]
        n_neg_test = neg_test_vectors.shape[0]
        logging.info(
            "Setting up RF classification with pos train (difference vectors): {}, neg train {}, pos test {}, neg test {}"
                .format(n_pos_train, n_neg_train, n_pos_test, n_neg_test))
        X_test = pd.concat([pos_test_vectors, neg_test_vectors])
        y_test = np.concatenate((np.ones(n_pos_test), np.zeros(n_neg_test)))
        # Perform random grid search for best parameters using the training data
        random_grid = KcetRandomForest._init_random_grid()
        rf = RandomForestClassifier()
        rf_random = RandomizedSearchCV(estimator=rf, param_distributions=random_grid, n_iter=1, cv=10, random_state=42)
        rf_random.fit(X_train, y_train)
        best_model = rf_random.best_estimator_
        # Now estimate the performance on the held out testing data
        y_pred = best_model.predict(X_test)
        yproba = best_model.predict_proba(X_test)[::, 1]
        return y_pred, y_test, yproba, n_pos_train, n_neg_train, n_pos_test, n_neg_test

    @staticmethod
    def _init_random_grid():
        """
        This is a convenience method that specifies the range of parameters we will
        use for the random grid search to optimize the random forest
        """
        # Number of trees in random forest
        n_estimators = [100, 200, 300, 400, 500]
        # Number of features to consider at every split
        max_features = ['auto', 'sqrt']
        # Maximum number of levels in tree
        max_depth = [10, 20, 30, 40, 50, None]
        # Minimum number of samples required to split a node
        min_samples_split = [2, 3, 5, 7, 10]
        # Minimum number of samples required at each leaf node
        min_samples_leaf = [1, 2, 4]
        # Method of selecting samples for training each tree
        bootstrap = [True, False]
        # Create the random grid
        random_grid = {'n_estimators': n_estimators,
                       'max_features': max_features,
                       'max_depth': max_depth,
                       'min_samples_split': min_samples_split,
                       'min_samples_leaf': min_samples_leaf,
                       'bootstrap': bootstrap}
        return random_grid
