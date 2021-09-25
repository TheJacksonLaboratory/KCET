from kcet.kcet_dataset_generator import KcetDatasetGenerator
import os
from unittest import TestCase


class TestKCETDatasetGenerator(TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Set up small fake dataset for testing
        """
        current_dir = os.path.dirname(__file__)
        ct_by_phase_path_1 = os.path.join(current_dir, 'data', 'small_ct_by_phase.tsv')
        embeddings = os.path.join(current_dir, 'data', 'embeddings10.npy')
        words = os.path.join(current_dir, 'data', 'words10.txt')
        cls.data_generator = KcetDatasetGenerator(clinical_trials=ct_by_phase_path_1, embeddings=embeddings, words=words)

    def test_get_positive_training_data_set_1_2015(self):
        """
        Note that we only regard phase 4 as 'positive' for training purposes
        There is one drug-disease link phase 4 up to 2015:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174

        afatinib targets 3 protein kinases:
        afatinib	EGFR
        afatinib	ERBB2
        afatinib    ERBB4

        So, there are 3 links in positive validation set:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB4   ncbigene2066
        """
        df_pos_training = self.data_generator._get_positive_training_data_set(2015)
        self.assertEqual(3, df_pos_training.shape[0])

    def test_get_positive_training_data_2006(self):
        """
        There is no drug-disease link phase 4 up to 2006:
        """
        df_pos_training = self.data_generator._get_positive_training_data_set(2006)
        self.assertEqual(0, df_pos_training.shape[0])



    def test_get_positive_validation_data_set_1_phase4_later_year_2012_2015(self):
        """
        There is one disease-drug link between 2012 and 2015:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174

        But,
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        was studied in phase 1,2,3 clinical trials studies.



        So, there is no  kinase-cancer link

        """
        df_pos_validation = self.data_generator.get_positive_test_embeddings(2012, 2013, 2015, phase4=True)
        #print(df_pos_validation)
        self.assertEqual(0, df_pos_validation.shape[0])

    def test_get_positive_validation_data_set_1_phase4_2012_2014_2015(self):
        """
        There is one disease-drug link between 2012 and 2015:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174

        But,
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        was studied in phase 1,2,3 clinical trials studies.



        So, there is no  kinase-cancer link

        """
        df_pos_validation = self.data_generator.get_positive_test_embeddings(2012, 2014, 2015, phase4=True)
        self.assertEqual(0, df_pos_validation.shape[0])

    def test_get_positive_validation_data_set_1_phase4_2012_2014_2020(self):
        """
        There are is one disease-drug link between 2014 and 2020 of phase 4:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174


        But,
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        was studied in phase 1,2,3 clinical trials studies.

        """
        df_pos_validation = self.data_generator.get_positive_test_embeddings(2012, 2014, 2020, phase4=True)
        self.assertEqual(0, df_pos_validation.shape[0])

    def test_get_positive_validation_data_set_1_phase4_later_year_2013_2018(self):
        """
         There is one disease-drug link between 2012 and 2015:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174

        But,
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        was studied in phase 1,2,3 clinical trials studies.

        So, there is no  kinase-cancer link

        """
        df_pos_validation = self.data_generator.get_positive_test_embeddings(2013, 2014, 2018)
        # print(df_pos_validation)
        self.assertEqual(0, df_pos_validation.shape[0])

    def test_get_negative_training_dataset_1_2008(self):
        """
        There are 2 lines up to 2008.
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 2	2007	2020	NCT02716311;NCT02098954;NCT02369484;NCT03810872;NCT02595840;NCT02795156;NCT03574402;NCT02747953;NCT02488694;NCT03157089;NCT04148898;NCT03623750;NCT02470065;NCT02597946;NCT04470076;NCT01932229;NCT01542437;NCT03727724;NCT04497584;NCT01415011;NCT02906163;NCT01003899;NCT02183883;NCT00730925;NCT00525148;NCT03399669;NCT00711594;NCT01156545;NCT00796549;NCT02450656;NCT01746251
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 3	2008	2015	NCT02044380;NCT02438722;NCT01853826;NCT01523587;NCT01814553;NCT01121393;NCT00656136;NCT00949650;NCT01085136;NCT01953913

        But, there is no disease-drug of phase 4.


        Number of negative training links is 10 times the number of positive training set. So, there are 0 negative training links
        """
        df_pos_training = self.data_generator._get_positive_training_data_set(2008)
        df_neg_training = self.data_generator.get_neg_training_embeddings(target_year=2008, n_neg_examples=10*len(df_pos_training))
        self.assertEqual(0, df_neg_training.shape[0])
