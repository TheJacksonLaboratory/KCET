from kcet.ct_by_phase_parser import CTParserByPhase
import os
from unittest import TestCase


class TestParseCTbyPhase(TestCase):

    @classmethod
    def setUpClass(cls):
        current_dir = os.path.dirname(__file__)
        ct_by_phase_path = os.path.join(current_dir, 'data', "small_ct_by_phase.tsv")
        parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
        drug_kinase_links_data_path = os.path.join(parent_dir, 'input', 'drug_kinase_links.tsv')

        phase = 4
        date = 2020
        cls.ct_parser = CTParserByPhase(clinical_trials=ct_by_phase_path,
                                        drug_kinase_links_data_path=drug_kinase_links_data_path,
                                        phase=phase,
                                        date=date)
        cls.ct_parser.get_kinase_cancer_links()
        cls.ct_parser.num_kinase_cancer_links()

    def test_num_kinase_disease_links(self):
        """
        Test the number of kinase-disease links by parsing small_ct_by_phase.tsv and drug_kinase_links.tsv.
        There are 1 line of data with Phase 4 and start_date before 2020:
        Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 4	2014	2020	NCT04413201;NCT02695290;NCT02208843;NCT04356118;NCT02514174
        In drug_kinase_links.tsv, for the drug afatinib, there are 2 lines:
        afatinib	EGFR
        afatinib	ERBB2

        So, there are two kinase-disease_links:
        Carcinoma, Non-Small-Cell Lung EGFR
        Carcinoma, Non-Small-Cell Lung ERBB2
        """
        expected_num_kinase_disease_links = self.ct_parser.num_kinase_cancer_links()
        self.assertEqual(expected_num_kinase_disease_links, 2)
