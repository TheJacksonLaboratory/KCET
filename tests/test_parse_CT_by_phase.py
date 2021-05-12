from kcet.ct_by_phase_parser import CTParserByPhase
import os
from unittest import TestCase


class TestParseCTbyPhase(TestCase):

    @classmethod
    def setUpClass(cls):
        current_dir = os.path.dirname(__file__)
        ct_by_phase_path = os.path.join(current_dir, 'data', "small_ct_by_phase.tsv")
        target_year = 2020
        cls.ct_parser = CTParserByPhase(clinical_trials=ct_by_phase_path,
                                        year=target_year)

    def test_num_disease_cancer_links(self):
        """
        Test the number of kinase-disease links by parsing small_ct_by_phase.tsv.
        Urethral Neoplasms-afatinib
        Ureteral Neoplasms-afatinib
        Multiple Myeloma-afatinib
        Urinary Bladder Neoplasms-afatinib
        Carcinoma, Non-Small-Cell Lung-afatinib


        So, there are two kinase-disease_links:
        Carcinoma, Non-Small-Cell Lung EGFR
        Carcinoma, Non-Small-Cell Lung ERBB2
        """
        disease_links = self.ct_parser.get_all_phases()
        seen = set()
        for index, row in disease_links.iterrows():
            cancer = row['cancer']
            pki = row['pki']
            key = "%s-%s" % (cancer, pki)
            seen.add(key)
        self.assertEqual(5, len(seen))

    def test_num_disease_cancer_links_phase_4(self):
        """
        Test the number of kinase-disease links of phase 4 by parsing small_ct_by_phase.tsv.
        Urethral Neoplasms-afatinib
        Ureteral Neoplasms-afatinib
        Multiple Myeloma-afatinib
        Urinary Bladder Neoplasms-afatinib
        Carcinoma, Non-Small-Cell Lung-afatinib


        So, there is only one kinase-disease_link of pahse 4:
        arcinoma, Non-Small-Cell Lung - afatinib
        """
        disease_links = self.ct_parser.get_phase_4()
        seen = set()
        for index, row in disease_links.iterrows():
            cancer = row['cancer']
            pki = row['pki']
            key = "%s-%s" % (cancer, pki)
            seen.add(key)
        self.assertEqual(1, len(seen))