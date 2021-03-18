from kcet.kcet_dataset_generator import KcetDatasetGenerator
import os
from unittest import TestCase


class TestKCETDatasetGEnerator(TestCase):
    @classmethod
    def setUpClass(cls):
        current_dir = os.path.dirname(__file__)
        ct_by_phase_path_1 = os.path.join(current_dir, 'data', "small_ct_by_phase.tsv")
        cls.kcet_data_generator_1 = KcetDatasetGenerator(clinical_trials=ct_by_phase_path_1)

        ct_by_phase_path_2 = os.path.join(current_dir, 'data', "small_ct_by_phase_2.tsv")
        cls.kcet_data_generator_2 = KcetDatasetGenerator(clinical_trials=ct_by_phase_path_2)

    def test_get_positive_training_data_set_1_2015(self):
        """
        There is one drug-disease link phase 4 up to 2015:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174

        afatinib targets two protein kinases:
        afatinib	EGFR
        afatinib	ERBB2

        So, there are two links in positive validation set:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064

        """
        df_pos_training = self.kcet_data_generator_1._get_positive_data_set(2015)
        self.assertEqual(2, df_pos_training.shape[0])

    def test_get_positive_validation_data_set_1_2015(self):
        """
        There is one drug-disease links after 2015:
        Multiple Myeloma	D009101	afatinib	Phase 1	2020	2020	NCT03878524

        afatinib targets two protein kinases:
        afatinib	EGFR    ncbigene1956
        afatinib	ERBB2   ncbigene2064

        So, there are two links in positive validation set:
        Multiple Myeloma	D009101 EGFR    ncbigene1956
        Multiple Myeloma    D009101 ERBB2   ncbigene2064

        """
        pos_validation = self.kcet_data_generator_1._get_positive_validation_data_set(2015)
        #print(pos_validation)
        self.assertEqual(2, len(pos_validation))

    def test_get_positive_validation_data_set_1_2007(self):
        """
        There are 8 drug-disease links after 2007:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 1	2009	2020	NCT03711422;NCT03827070;NCT02191891;NCT01288430;NCT01647711;NCT03054038;NCT01090011;NCT00993499;NCT04448379;NCT01999985;NCT02364609
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 3	2008	2015	NCT02044380;NCT02438722;NCT01853826;NCT01523587;NCT01814553;NCT01121393;NCT00656136;NCT00949650;NCT01085136;NCT01953913
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        Urinary Bladder Neoplasms	D001749	afatinib	Phase 2	2013	2015	NCT02122172;NCT02465060
        Urethral Neoplasms	D014523	afatinib	Phase 2	2013	2013	NCT02122172
        Ureteral Neoplasms	D014516	afatinib	Phase 2	2013	2013	NCT02122172
        Multiple Myeloma	D009101	afatinib	Phase 1	2020	2020	NCT03878524
        Multiple Myeloma	D009101	afatinib	Phase 2	2015	2016	NCT02693535;NCT04439136;NCT02465060

        afatinib targets two protein kinases:
        afatinib	EGFR
        afatinib	ERBB2

        So, there are 10 links in positive validation set:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        Urinary Bladder Neoplasms	D001749 EGFR    ncbigene1956
        Urinary Bladder Neoplasms	D001749 ERBB2   ncbigene2064
        Urethral Neoplasms	D014523 EGFR    ncbigene1956
        Urethral Neoplasms	D014523 ERBB2   ncbigene2064
        Ureteral Neoplasms	D014516 EGFR    ncbigene1956
        Ureteral Neoplasms	D014516 ERBB2   ncbigene2064
        Multiple Myeloma	D009101 EGFR    ncbigene1956
        Multiple Myeloma    D009101 ERBB2   ncbigene2064

        """
        pos_validation = self.kcet_data_generator_1._get_positive_validation_data_set(2007)
        #print(pos_validation)
        self.assertEqual(10, pos_validation.shape[0])

    def test_get_positive_validation_data_set_1_phase_4_2015(self):
        """
        There is no drug-disease link of phase 4 after 2015
        So, there is no link in the positive validation set
        """
        df_pos_validation = self.kcet_data_generator_1._get_positive_validation_data_set_phase_4(2015)
        self.assertEqual(0, df_pos_validation.shape[0])

    def test_get_positive_validation_data_set_1_phase_4_2007(self):
        """
        There is 1 drug-disease link of phase 4 after 2007
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        afatinib targets two protein kinases:
        afatinib	EGFR
        afatinib	ERBB2
        So, there are two links in positive validation set:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        """
        df_pos_validation = self.kcet_data_generator_1._get_positive_validation_data_set_phase_4(2007)
        self.assertEqual(2, df_pos_validation.shape[0])

    def test_get_positive_validation_data_set_1_later_year_2014_2018(self):
        """
        There are three disease-drug links between 2014 and 2018:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        Multiple Myeloma	D009101	afatinib	Phase 1	2020	2020	NCT03878524
        Multiple Myeloma	D009101	afatinib	Phase 2	2015	2016	NCT02693535;NCT04439136;NCT02465060

        afatinib targets two protein kinases: EGFR and ERBB2

        So, there will be 4 kinase-cancer links:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        Multiple Myeloma	D009101  EGFR   ncbigene1956
        Multiple Myeloma  D009101 ERBB2 ncbigene2064
        """
        df_pos_validation = self.kcet_data_generator_1._get_positive_validation_data_set_later_year(2013, 2018)
        #print(df_pos_validation)
        self.assertEqual(4, df_pos_validation.shape[0])



    def test_get_negative_training_dataset_1_2008(self):
        """
        There are 2 lines up to 2008.
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 2	2007	2020	NCT02716311;NCT02098954;NCT02369484;NCT03810872;NCT02595840;NCT02795156;NCT03574402;NCT02747953;NCT02488694;NCT03157089;NCT04148898;NCT03623750;NCT02470065;NCT02597946;NCT04470076;NCT01932229;NCT01542437;NCT03727724;NCT04497584;NCT01415011;NCT02906163;NCT01003899;NCT02183883;NCT00730925;NCT00525148;NCT03399669;NCT00711594;NCT01156545;NCT00796549;NCT02450656;NCT01746251
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 3	2008	2015	NCT02044380;NCT02438722;NCT01853826;NCT01523587;NCT01814553;NCT01121393;NCT00656136;NCT00949650;NCT01085136;NCT01953913


        There are two positive training link:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064

        Number of negative training links is 10 times the number of positive training set. So, there are 20 negative training links
        """
        df_neg_training = self.kcet_data_generator_1._get_negative_training_dataset(2008)
        self.assertEqual(20, df_neg_training.shape[0])

    def test_get_negative_validation_dataset_1_2008(self):
        """
        There are 2 lines up to 2008.
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 2	2007	2020	NCT02716311;NCT02098954;NCT02369484;NCT03810872;NCT02595840;NCT02795156;NCT03574402;NCT02747953;NCT02488694;NCT03157089;NCT04148898;NCT03623750;NCT02470065;NCT02597946;NCT04470076;NCT01932229;NCT01542437;NCT03727724;NCT04497584;NCT01415011;NCT02906163;NCT01003899;NCT02183883;NCT00730925;NCT00525148;NCT03399669;NCT00711594;NCT01156545;NCT00796549;NCT02450656;NCT01746251
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 3	2008	2015	NCT02044380;NCT02438722;NCT01853826;NCT01523587;NCT01814553;NCT01121393;NCT00656136;NCT00949650;NCT01085136;NCT01953913


        There are two positive training links:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064

        Number of negative training links is 10 times the number of positive training set. So, there are 20 negative training links
        We generate 20 negative validation set. The negative links should not exist
        """
        df_neg_training = self.kcet_data_generator_1._get_negative_training_dataset(2008)
        df_neg_validation = self.kcet_data_generator_1._get_negative_validation_data_set(df_neg_training,2008)
        #print(df_neg_validation)
        self.assertEqual(20, df_neg_validation.shape[0])

        ##################################################################################

    def test_get_positive_training_data_set_2_2015(self):
        """
        There is one drug-disease link phase 4 up to 2015:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174

        afatinib targets two protein kinases:
        afatinib	EGFR
        afatinib	ERBB2

        So, there are two links in positive validation set:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064

        """
        df_pos_training = self.kcet_data_generator_2._get_positive_data_set(2015)
        self.assertEqual(2, df_pos_training.shape[0])

    def test_get_positive_training_data_set_2_2020(self):
        """
        There are two drug-disease link phase 4 up to 2020:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        Breast Neoplasms	D001943	abemaciclib	Phase 4	2019	2021	NCT04707196;NCT03988114;NCT04031885

        afatinib targets two protein kinases:
        afatinib	EGFR
        afatinib	ERBB2
        abemaciclib targets two protein kinases:
        abemaciclib CDK4
        abemaciclib CDK6

        So, there are 4 links in positive validation set:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        Breast Neoplasms	D001943 CDK4 ncbigene1019
        Breast Neoplasms	D001943 CDK6 ncbigene1021
        """
        df_pos_training = self.kcet_data_generator_2._get_positive_data_set(2020)
        self.assertEqual(4, df_pos_training.shape[0])