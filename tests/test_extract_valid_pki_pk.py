from kcet import DrugCentralPkPkiParser
from unittest import TestCase
import logging
logging.disable(logging.CRITICAL)

class TestExtractValidPkiPk(TestCase):
    @classmethod
    def setUpClass(cls):
        pkpki = DrugCentralPkPkiParser()
        cls.valid_pki_pk = pkpki.get_pk_pki_with_threshold()
        cls.EPSILON = 0.00001

    def test_afatinib_EGFR(self):
        """
        afatinib	EGFR	18408761	0.0001	Kd	INHIBITOR
        """
        df = self.valid_pki_pk
        item = df.loc[(df['PKI'] == 'afatinib') & (df['PK'] == 'EGFR')]
        self.assertEqual(1, len(item))
        row = item.iloc[0]
        kd = float(row['ACT_VALUE'])
        pmid = row['PMID']
        self.assertAlmostEqual(kd, 0.0001)
        self.assertEqual(pmid, '18408761')

    def test_alpelisib_PIK3CA(self):
        """
        alpelisib	PIK3CA	25544637
        alpelisib	PIK3CA	n/a	0.004600446663935	IC50
        """
        df = self.valid_pki_pk
        item = df.loc[(df['PKI'] == 'alpelisib') & (df['PK'] == 'PIK3CA')]
        self.assertEqual(2, len(item))
        row = item.iloc[0]
        medication = row['PKI']
        kd = row['ACT_VALUE']
        pmid = row['PMID']
        self.assertEqual('alpelisib', medication)
        self.assertEqual(kd, 'n/a')  # We do not have Kd data for this combination
        self.assertEqual(pmid, '25544637')

    def test_apatinib_KIT_removed(self):
        """
        This entry should have been removed since the Kd is too high
        apatinib	KIT	20923544	0.428999964569677	IC50
        """
        df = self.valid_pki_pk
        items = df.loc[(df['PKI'] == 'apatinib') & (df['PK'] == 'KIT')]
        self.assertEqual(0, len(items.index))
