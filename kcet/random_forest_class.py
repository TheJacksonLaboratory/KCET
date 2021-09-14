from .kcet_dataset_generator import KcetDatasetGenerator
from .kinase_predictor import KinasePredictor

import pandas as pd
import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV

class KcetRandomForest:
    """
    This class is a wrapper around scikit learn functions for random forest classification.
    """

    def __init__(self, 
                data_gen: KcetDatasetGenerator,
                embedddingfile: str,
                wordsfile: str,
                target: int) -> None:

        self._data_generator = data_gen
        self._target_year = target
        if not os.path.isfile(embedddingfile):
            raise FileNotFoundError("Could not find embedding file at " + embedddingfile)
        if not os.path.isfile(wordsfile):
            raise FileNotFoundError("Could not find embedding/words file at " + wordsfile)
        self._predictor = KinasePredictor(embeddings=embedddingfile, words=wordsfile)
    

    def classify(self, num_years_later: int, mid_year: int):
        """
        Perform random forest learning. From the vectors extracted from the data from the target year, predict
        clinical trials starting at midyear and going num_years_later
        For instance, 
        target_year = 2010
        mid_year = 2019
        num_years_after_the_mid_year = 1
        creates test tests from 2019 to 2020.
        """
        pos_train_df, neg_train_df, pos_validation_df, neg_validation_df = \
            self._data_generator.get_data_years_after_target_year_upto_later_year(target_year=self._target_year, mid_year = mid_year, num_years_later= num_years_later)
        #diff_vectors_pos_validation = pd.read_pickle(positive_validation_pickle_path)
        #diff_vectors_neg_validation = pd.read_pickle(negative_validation_pickle_path)
        #diff_vectors_pos_training = pd.read_pickle(positive_train_pickle_path)
        #diff_vectors_neg_training = pd.read_pickle(negative_train_pickle_path)
        # The following gets the names of the vectors
        X_train = pd.concat([pos_train_df, neg_train_df])
        n_pos_train = pos_train_df.shape[0]
        n_neg_train = neg_train_df.shape[0]
        print(pos_train_df.head())

        #print("Total training vectors: %d" % len(X_train))
        y_train = np.concatenate((np.ones(n_pos_train), np.zeros(n_neg_train)))
        X_test = pd.concat([pos_validation_df, neg_validation_df])
        # print("Total test vectors: %d" % len(X_test))
        n_pos_test = pos_validation_df.shape[0]
        n_neg_test = neg_validation_df.shape[0]
        y_test = np.concatenate((np.ones(n_pos_test), np.zeros(n_neg_test)))
        # print("Total test labels: %d" % len(y_test))
        random_grid = KcetRandomForest._init_random_grid()

        rf = RandomForestClassifier()
        rf_random = RandomizedSearchCV(estimator = rf, param_distributions = random_grid, n_iter = 1, cv = 10, random_state=42)
        rf_random.fit(X_train,y_train)
        best_model = rf_random.best_estimator_

        y_pred = best_model.predict(X_test)
        yproba = best_model.predict_proba(X_test)[::,1]


    def _get_positive_example_vectors(self):



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
        #Maximum number of levels in tree
        max_depth = [10, 20, 30, 40, 50]
        max_depth.append(None)
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