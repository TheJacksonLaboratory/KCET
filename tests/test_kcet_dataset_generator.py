from kcet.kcet_dataset_generator import KcetDatasetGenerator
import os
from unittest import TestCase


class TestKCETDatasetGEnerator(TestCase):
    @classmethod
    def setUpClass(cls):
        current_dir = os.path.dirname(__file__)
        ct_by_phase_path = os.path.join(current_dir, 'data', "small_ct_by_phase.tsv")
        cls.kcet_data_generator = KcetDatasetGenerator(clinical_trials=ct_by_phase_path)

    def test_get_positive_validation_data_set(self):
        """
        There is one drug-disease links after 2015:
        Multiple Myeloma	D009101	afatinib	Phase 1	2020	2020	NCT03878524

        afatinib targets two protein kinases:
        afatinib	EGFR	18408761
        afatinib	ERBB2	18408761

        So, there are two links in positive validation set:
        Multiple Myeloma	D009101 EGFR
        Multiple Myeloma    D009101 ERBB2

        """
        df_pos_validation = self.kcet_data_generator._get_positive_validation_data_set(2015)
        self.assertEqual(2, df_pos_validation.shape[0])

    def test_get_positive_validation_data_set_phase_4(self):
        """
        There is no drug-disease link of phase 4 after 2015
        So, there is no link in the positive validation set
        """
        df_pos_validation = self.kcet_data_generator._get_positive_validation_data_set_phase_4(2015)
        self.assertEqual(0, df_pos_validation.shape[0])

    def test_get_positive_validation_data_set_later_year(self):
        """
        There are three disease-drug links between 2014 and 2018:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        Multiple Myeloma	D009101	afatinib	Phase 1	2020	2020	NCT03878524
        Multiple Myeloma	D009101	afatinib	Phase 2	2015	2016	NCT02693535;NCT04439136;NCT02465060

        afatinib targets two protein kinases: EGFR and ERBB2

        So, there will be 4 kinase-cancer links:
        Carcinoma, Non-Small-Cell Lung	D002289 EGFR
        Carcinoma, Non-Small-Cell Lung	D002289 ERBB2
        Multiple Myeloma	D009101  EGFR
        Multiple Myeloma  D009101 ERBB2
        """
        df_pos_validation = self.kcet_data_generator._get_positive_validation_data_set_later_year(2013, 2018)
        self.assertEqual(4, df_pos_validation.shape[0])