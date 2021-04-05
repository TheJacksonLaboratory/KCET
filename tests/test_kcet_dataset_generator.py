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

        afatinib targets 3 protein kinases:
        afatinib	EGFR
        afatinib	ERBB2
        afatinib    ERBB4

        So, there are 3 links in positive validation set:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB4   ncbigene2066


        """
        df_pos_training = self.kcet_data_generator_1._get_positive_training_data_set(2015)
        self.assertEqual(3, df_pos_training.shape[0])

    def test_get_positive_training_data_set_1_2006(self):
        """
        There is no drug-disease link phase 4 up to 2006:
        """
        df_pos_training = self.kcet_data_generator_1._get_positive_training_data_set(2006)
        print(df_pos_training.head())
        self.assertEqual(0, df_pos_training.shape[0])

    def test_get_positive_validation_data_set_1_later_year_2013_2018(self):
        """
        There are three disease-drug links between 2014 and 2018:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        Multiple Myeloma	D009101	afatinib	Phase 1	2020	2020	NCT03878524
        Multiple Myeloma	D009101	afatinib	Phase 2	2015	2016	NCT02693535;NCT04439136;NCT02465060

        afatinib targets 3 protein kinases: EGFR, ERBB2, ERBB4

        So, there will be 6 kinase-cancer links:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB4  ncbigene2066
        Multiple Myeloma	D009101  EGFR   ncbigene1956
        Multiple Myeloma  D009101 ERBB2 ncbigene2064
        Multiple Myeloma  D009101 ERBB4 ncbigene2066
        """
        df_pos_validation = self.kcet_data_generator_1._get_positive_validation_data_set_later_year(2013, 2018)
        #print(df_pos_validation)
        self.assertEqual(6, df_pos_validation.shape[0])

    def test_get_positive_validation_data_set_1_phase4_later_year_2013_2018(self):
        """
        There is only one disease-drug link between 2014 and 2018:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174

        afatinib targets 3 protein kinases: EGFR, ERBB2, ERBB4

        So, there will be 6 kinase-cancer links:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB4  ncbigene2066
        """
        df_pos_validation = self.kcet_data_generator_1._get_positive_validation_data_set_later_year_phase_4(2013, 2018)
        # print(df_pos_validation)
        self.assertEqual(3, df_pos_validation.shape[0])

    def test_get_negative_training_dataset_1_2008(self):
        """
        There are 2 lines up to 2008.
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 2	2007	2020	NCT02716311;NCT02098954;NCT02369484;NCT03810872;NCT02595840;NCT02795156;NCT03574402;NCT02747953;NCT02488694;NCT03157089;NCT04148898;NCT03623750;NCT02470065;NCT02597946;NCT04470076;NCT01932229;NCT01542437;NCT03727724;NCT04497584;NCT01415011;NCT02906163;NCT01003899;NCT02183883;NCT00730925;NCT00525148;NCT03399669;NCT00711594;NCT01156545;NCT00796549;NCT02450656;NCT01746251
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 3	2008	2015	NCT02044380;NCT02438722;NCT01853826;NCT01523587;NCT01814553;NCT01121393;NCT00656136;NCT00949650;NCT01085136;NCT01953913

        But, there is no disease-drug of phase 4.


        Number of negative training links is 10 times the number of positive training set. So, there are 0 negative training links
        """
        df_pos_training = self.kcet_data_generator_1._get_positive_training_data_set(2008)
        df_neg_training = self.kcet_data_generator_1._get_negative_training_dataset(df_pos_training,2008)
        self.assertEqual(0, df_neg_training.shape[0])

    def test_get_negative_validation_dataset_1_2014(self):
        """
        There are 7 lines up to 2014.
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 1	2009	2020	NCT03711422;NCT03827070;NCT02191891;NCT01288430;NCT01647711;NCT03054038;NCT01090011;NCT00993499;NCT04448379;NCT01999985;NCT02364609
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 2	2007	2020	NCT02716311;NCT02098954;NCT02369484;NCT03810872;NCT02595840;NCT02795156;NCT03574402;NCT02747953;NCT02488694;NCT03157089;NCT04148898;NCT03623750;NCT02470065;NCT02597946;NCT04470076;NCT01932229;NCT01542437;NCT03727724;NCT04497584;NCT01415011;NCT02906163;NCT01003899;NCT02183883;NCT00730925;NCT00525148;NCT03399669;NCT00711594;NCT01156545;NCT00796549;NCT02450656;NCT01746251
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 3	2008	2015	NCT02044380;NCT02438722;NCT01853826;NCT01523587;NCT01814553;NCT01121393;NCT00656136;NCT00949650;NCT01085136;NCT01953913
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        Urinary Bladder Neoplasms	D001749	afatinib	Phase 2	2013	2015	NCT02122172;NCT02465060
        Urethral Neoplasms	D014523	afatinib	Phase 2	2013	2013	NCT02122172
        Ureteral Neoplasms	D014516	afatinib	Phase 2	2013	2013	NCT02122172

        There are only 3 positive links of phase 4 up to 2014:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB4   ncbigene2066


        Number of negative training links is 10 times the number of positive training set. So, there are 30 negative training links
        We generate 30 negative validation set.
        """
        df_pos_training = self.kcet_data_generator_1._get_positive_training_data_set(2014)
        df_neg_training = self.kcet_data_generator_1._get_negative_training_dataset(df_pos_training,2014)
        df_neg_validation = self.kcet_data_generator_1._get_negative_validation_data_set(df_neg_training,2014)
        #print(df_neg_validation)
        self.assertEqual(30, df_neg_validation.shape[0])

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
        df_pos_training = self.kcet_data_generator_2._get_positive_training_data_set(2015)
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
        df_pos_training = self.kcet_data_generator_2._get_positive_training_data_set(2020)
        self.assertEqual(4, df_pos_training.shape[0])

    def test_get_positive_validation_data_set_2_2015(self):
        """
        There are two drug-disease links after 2015:
        Multiple Myeloma	D009101	afatinib	Phase 1	2020	2020	NCT03878524
        Breast Neoplasms	D001943	abemaciclib	Phase 4	2019	2021	NCT04707196;NCT03988114;NCT04031885

        afatinib targets two protein kinases:
        afatinib	EGFR    ncbigene1956
        afatinib	ERBB2   ncbigene2064

        abemaciclib targets two protein kinases:
        abemaciclib CDK4 ncbigene1019
        abemaciclib CDK6 ncbigene1021

        So, there are 4 links in positive validation set:
        Multiple Myeloma	D009101 EGFR    ncbigene1956
        Multiple Myeloma    D009101 ERBB2   ncbigene2064
        Breast Neoplasms	D001943 CDK4    ncbigene1019
        Breast Neoplasms	D001943 CDK6    ncbigene1021

        """
        pos_validation = self.kcet_data_generator_2._get_positive_validation_data_set(2015)
        print(pos_validation)
        self.assertEqual(4, len(pos_validation))

    def test_get_positive_validation_data_set_2_2007(self):
        """
        There are 12 drug-disease links after 2007:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 1	2009	2020	NCT03711422;NCT03827070;NCT02191891;NCT01288430;NCT01647711;NCT03054038;NCT01090011;NCT00993499;NCT04448379;NCT01999985;NCT02364609
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 3	2008	2015	NCT02044380;NCT02438722;NCT01853826;NCT01523587;NCT01814553;NCT01121393;NCT00656136;NCT00949650;NCT01085136;NCT01953913
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        Urinary Bladder Neoplasms	D001749	afatinib	Phase 2	2013	2015	NCT02122172;NCT02465060
        Urethral Neoplasms	D014523	afatinib	Phase 2	2013	2013	NCT02122172
        Ureteral Neoplasms	D014516	afatinib	Phase 2	2013	2013	NCT02122172
        Multiple Myeloma	D009101	afatinib	Phase 1	2020	2020	NCT03878524
        Multiple Myeloma	D009101	afatinib	Phase 2	2015	2016	NCT02693535;NCT04439136;NCT02465060
        Breast Neoplasms	D001943	abemaciclib	Phase 1	2014	2021	NCT04188548;NCT04316169;NCT04481113;NCT04483505;NCT02779751;NCT03099174;NCT04585724;NCT04514159;NCT03846583;NCT03616587;NCT02057133;NCT03878524;NCT04088032
        Breast Neoplasms	D001943	abemaciclib	Phase 2	2014	2021	NCT04523857;NCT04352777;NCT04305834;NCT02747004;NCT04432454;NCT04256941;NCT03979508;NCT04603183;NCT04351230;NCT02308020;NCT03939897;NCT03130439;NCT03227328;NCT02831530;NCT03280563;NCT02675231;NCT03703466;NCT04614194;NCT04293393;NCT02102490;NCT04305236;NCT02441946;NCT03913234;NCT04227327
        Breast Neoplasms	D001943	abemaciclib	Phase 3	2014	2021	NCT02763566;NCT03155997;NCT04565054;NCT04752332;NCT03425838;NCT02246621;NCT04158362;NCT02107703
        Breast Neoplasms	D001943	abemaciclib	Phase 4	2019	2021	NCT04707196;NCT03988114;NCT04031885

        afatinib targets two protein kinases:
        afatinib	EGFR    ncbigene1956
        afatinib	ERBB2   ncbigene2064

        abemaciclib targets two protein kinases:
        abemaciclib CDK4 ncbigene1019
        abemaciclib CDK6 ncbigene1021

        So, there are 12 links in positive validation set:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        Urinary Bladder Neoplasms	D001749 EGFR    ncbigene1956
        Urinary Bladder Neoplasms	D001749 ERBB2   ncbigene2064
        Urethral Neoplasms	D014523 EGFR    ncbigene1956
        Urethral Neoplasms	D014523 ERBB2   ncbigene2064
        Ureteral Neoplasms	D014516 EGFR    ncbigene1956
        Ureteral Neoplasms	D014516 ERBB2   ncbigene2064
        Multiple Myeloma	D009101 EGFR    ncbigene1956
        Multiple Myeloma	D009101 ERBB2   ncbigene2064
        Breast Neoplasms	D001943 CDK4    ncbigene1019
        Breast Neoplasms	D001943 CDK6    ncbigene1021

        """
        pos_validation = self.kcet_data_generator_2._get_positive_validation_data_set(2007)
        #print(pos_validation)
        self.assertEqual(12, len(pos_validation))

    def test_get_positive_validation_data_set_phase4_2_2015(self):
        """
        There is one drug-disease links of phase 4 after 2015:
        Breast Neoplasms	D001943	abemaciclib	Phase 4	2019	2021	NCT04707196;NCT03988114;NCT04031885



        abemaciclib targets two protein kinases:
        abemaciclib CDK4 ncbigene1019
        abemaciclib CDK6 ncbigene1021

        So, there are 2 links in positive validation set:
        Breast Neoplasms	D001943 CDK4    ncbigene1019
        Breast Neoplasms	D001943 CDK6    ncbigene1021

        """
        pos_validation = self.kcet_data_generator_2._get_positive_validation_data_set_phase_4(2015)
        print(pos_validation)
        self.assertEqual(2, len(pos_validation))

    def test_get_positive_validation_data_set_phase4_2_2007(self):
        """
        There are 12 drug-disease links after 2007:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        Breast Neoplasms	D001943	abemaciclib	Phase 4	2019	2021	NCT04707196;NCT03988114;NCT04031885

        afatinib targets two protein kinases:
        afatinib	EGFR    ncbigene1956
        afatinib	ERBB2   ncbigene2064

        abemaciclib targets two protein kinases:
        abemaciclib CDK4 ncbigene1019
        abemaciclib CDK6 ncbigene1021

        So, there are 4 links in positive validation set:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        Breast Neoplasms	D001943 CDK4    ncbigene1019
        Breast Neoplasms	D001943 CDK6    ncbigene1021

        """
        pos_validation = self.kcet_data_generator_2._get_positive_validation_data_set_phase_4(2007)
        #print(pos_validation)
        self.assertEqual(4, len(pos_validation))

    def test_get_positive_validation_data_set_2_later_year_2013_2018(self):
        """
        There are 6 disease-drug links between 2014 and 2018:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        Multiple Myeloma	D009101	afatinib	Phase 1	2020	2020	NCT03878524
        Multiple Myeloma	D009101	afatinib	Phase 2	2015	2016	NCT02693535;NCT04439136;NCT02465060
        Breast Neoplasms	D001943	abemaciclib	Phase 1	2014	2021	NCT04188548;NCT04316169;NCT04481113;NCT04483505;NCT02779751;NCT03099174;NCT04585724;NCT04514159;NCT03846583;NCT03616587;NCT02057133;NCT03878524;NCT04088032
        Breast Neoplasms	D001943	abemaciclib	Phase 2	2014	2021	NCT04523857;NCT04352777;NCT04305834;NCT02747004;NCT04432454;NCT04256941;NCT03979508;NCT04603183;NCT04351230;NCT02308020;NCT03939897;NCT03130439;NCT03227328;NCT02831530;NCT03280563;NCT02675231;NCT03703466;NCT04614194;NCT04293393;NCT02102490;NCT04305236;NCT02441946;NCT03913234;NCT04227327
        Breast Neoplasms	D001943	abemaciclib	Phase 3	2014	2021	NCT02763566;NCT03155997;NCT04565054;NCT04752332;NCT03425838;NCT02246621;NCT04158362;NCT02107703

        afatinib targets two protein kinases:
        afatinib	EGFR    ncbigene1956
        afatinib	ERBB2   ncbigene2064

        abemaciclib targets two protein kinases:
        abemaciclib CDK4 ncbigene1019
        abemaciclib CDK6 ncbigene1021

        So, there will be 6 kinase-cancer links:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        Multiple Myeloma	D009101  EGFR   ncbigene1956
        Multiple Myeloma  D009101 ERBB2 ncbigene2064
        Breast Neoplasms	D001943 CDK4    ncbigene1019
        Breast Neoplasms	D001943 CDK6    ncbigene1021
        """
        df_pos_validation = self.kcet_data_generator_2._get_positive_validation_data_set_later_year(2013, 2018)
        #print(df_pos_validation)
        self.assertEqual(6, df_pos_validation.shape[0])

    def test_get_positive_validation_data_set_2_later_year_2013_2020(self):
        """
        There are 6 disease-drug links between 2014 and 2020:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        Multiple Myeloma	D009101	afatinib	Phase 1	2020	2020	NCT03878524
        Multiple Myeloma	D009101	afatinib	Phase 2	2015	2016	NCT02693535;NCT04439136;NCT02465060
        Breast Neoplasms	D001943	abemaciclib	Phase 1	2014	2021	NCT04188548;NCT04316169;NCT04481113;NCT04483505;NCT02779751;NCT03099174;NCT04585724;NCT04514159;NCT03846583;NCT03616587;NCT02057133;NCT03878524;NCT04088032
        Breast Neoplasms	D001943	abemaciclib	Phase 2	2014	2021	NCT04523857;NCT04352777;NCT04305834;NCT02747004;NCT04432454;NCT04256941;NCT03979508;NCT04603183;NCT04351230;NCT02308020;NCT03939897;NCT03130439;NCT03227328;NCT02831530;NCT03280563;NCT02675231;NCT03703466;NCT04614194;NCT04293393;NCT02102490;NCT04305236;NCT02441946;NCT03913234;NCT04227327
        Breast Neoplasms	D001943	abemaciclib	Phase 3	2014	2021	NCT02763566;NCT03155997;NCT04565054;NCT04752332;NCT03425838;NCT02246621;NCT04158362;NCT02107703
        Breast Neoplasms	D001943	abemaciclib	Phase 4	2019	2021	NCT04707196;NCT03988114;NCT04031885

        afatinib targets two protein kinases:
        afatinib	EGFR    ncbigene1956
        afatinib	ERBB2   ncbigene2064

        abemaciclib targets two protein kinases:
        abemaciclib CDK4 ncbigene1019
        abemaciclib CDK6 ncbigene1021

        So, there will be 6 kinase-cancer links:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064
        Multiple Myeloma	D009101  EGFR   ncbigene1956
        Multiple Myeloma  D009101 ERBB2 ncbigene2064
        Breast Neoplasms	D001943 CDK4    ncbigene1019
        Breast Neoplasms	D001943 CDK6    ncbigene1021
        """
        df_pos_validation = self.kcet_data_generator_2._get_positive_validation_data_set_later_year(2013, 2020)
        #print(df_pos_validation)
        self.assertEqual(6, df_pos_validation.shape[0])

    def test_get_negative_training_dataset_2_2008(self):
        """
        There are 2 lines up to 2008.
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 2	2007	2020	NCT02716311;NCT02098954;NCT02369484;NCT03810872;NCT02595840;NCT02795156;NCT03574402;NCT02747953;NCT02488694;NCT03157089;NCT04148898;NCT03623750;NCT02470065;NCT02597946;NCT04470076;NCT01932229;NCT01542437;NCT03727724;NCT04497584;NCT01415011;NCT02906163;NCT01003899;NCT02183883;NCT00730925;NCT00525148;NCT03399669;NCT00711594;NCT01156545;NCT00796549;NCT02450656;NCT01746251
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 3	2008	2015	NCT02044380;NCT02438722;NCT01853826;NCT01523587;NCT01814553;NCT01121393;NCT00656136;NCT00949650;NCT01085136;NCT01953913

        But, there is no disease-drug link of phase 4 up to 2008.
        So, the number of positive links is 0.

        Number of negative training links is 10 times the number of positive training set. So, there are 0 negative training links
        """
        df_pos_training = self.kcet_data_generator_2._get_positive_training_data_set(2008)
        df_neg_training = self.kcet_data_generator_2._get_negative_training_dataset(df_pos_training,2008)
        self.assertEqual(0, df_neg_training.shape[0])

    def test_get_negative_training_dataset_2_2015(self):
        """
        There are 11 drug-disease links up to 2015:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 1	2009	2020	NCT03711422;NCT03827070;NCT02191891;NCT01288430;NCT01647711;NCT03054038;NCT01090011;NCT00993499;NCT04448379;NCT01999985;NCT02364609
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 2	2007	2020	NCT02716311;NCT02098954;NCT02369484;NCT03810872;NCT02595840;NCT02795156;NCT03574402;NCT02747953;NCT02488694;NCT03157089;NCT04148898;NCT03623750;NCT02470065;NCT02597946;NCT04470076;NCT01932229;NCT01542437;NCT03727724;NCT04497584;NCT01415011;NCT02906163;NCT01003899;NCT02183883;NCT00730925;NCT00525148;NCT03399669;NCT00711594;NCT01156545;NCT00796549;NCT02450656;NCT01746251
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 3	2008	2015	NCT02044380;NCT02438722;NCT01853826;NCT01523587;NCT01814553;NCT01121393;NCT00656136;NCT00949650;NCT01085136;NCT01953913
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        Urinary Bladder Neoplasms	D001749	afatinib	Phase 2	2013	2015	NCT02122172;NCT02465060
        Urethral Neoplasms	D014523	afatinib	Phase 2	2013	2013	NCT02122172
        Ureteral Neoplasms	D014516	afatinib	Phase 2	2013	2013	NCT02122172
        Multiple Myeloma	D009101	afatinib	Phase 2	2015	2016	NCT02693535;NCT04439136;NCT02465060
        Breast Neoplasms	D001943	abemaciclib	Phase 1	2014	2021	NCT04188548;NCT04316169;NCT04481113;NCT04483505;NCT02779751;NCT03099174;NCT04585724;NCT04514159;NCT03846583;NCT03616587;NCT02057133;NCT03878524;NCT04088032
        Breast Neoplasms	D001943	abemaciclib	Phase 2	2014	2021	NCT04523857;NCT04352777;NCT04305834;NCT02747004;NCT04432454;NCT04256941;NCT03979508;NCT04603183;NCT04351230;NCT02308020;NCT03939897;NCT03130439;NCT03227328;NCT02831530;NCT03280563;NCT02675231;NCT03703466;NCT04614194;NCT04293393;NCT02102490;NCT04305236;NCT02441946;NCT03913234;NCT04227327
        Breast Neoplasms	D001943	abemaciclib	Phase 3	2014	2021	NCT02763566;NCT03155997;NCT04565054;NCT04752332;NCT03425838;NCT02246621;NCT04158362;NCT02107703

        There is one drug-disease link of phase 4 up to 2015:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174

        afatinib targets two protein kinases:
        afatinib	EGFR    ncbigene1956
        afatinib	ERBB2   ncbigene2064


        There are 2 positive ilnks.
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR    ncbigene1956
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2   ncbigene2064

        Number of negative training links is 10 times the number of positive training set. So, there are 20 negative training links
        """
        df_pos_training = self.kcet_data_generator_2._get_positive_training_data_set(2015)
        df_neg_training = self.kcet_data_generator_2._get_negative_training_dataset(df_pos_training,2015)
        self.assertEqual(20, df_neg_training.shape[0])