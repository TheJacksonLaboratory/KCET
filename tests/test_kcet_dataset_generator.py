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
        df_pos_validation = self.kcet_data_generator._get_positive_validation_data_set(2015)
        self.assertEqual(2, df_pos_validation.shape[0])



    def test_get_positive_validation_data_set_phase_4(self):
        df_pos_validation = self.kcet_data_generator._get_positive_validation_data_set_phase_4(2015)
        self.assertEqual(0, df_pos_validation.shape[0])


    def test_get_positive_validation_data_set_later_year(self):
        df_pos_validation = self.kcet_data_generator._get_positive_validation_data_set_later_year(2015, 2018)
        self.assertEqual(0, df_pos_validation.shape[0])