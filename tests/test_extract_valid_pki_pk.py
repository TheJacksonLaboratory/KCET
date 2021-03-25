from kcet import PkPkiFilter
from unittest import TestCase
import logging
logging.disable(logging.CRITICAL)

class TestExtractValidPkiPk(TestCase):
    @classmethod
    def setUpClass(cls):
        pkpki = PkPkiFilter()
        cls.valid_pki_pk = pkpki.get_valid_pk_pki()
        cls.EPSILON = 0.00001

    def test_afatinib_EGFR(self):
        """
        afatinib	EGFR	18408761	0.0001	Kd	INHIBITOR
        """
        df = self.valid_pki_pk
        item = df.loc[(df['PKI'] == 'afatinib') & (df['PK'] == 'EGFR')]
        self.assertEqual(1, len(item))
        row = item.to_numpy()[0]
        self.assertAlmostEqual(0.0001, float(row[2]))
        self.assertEqual('18408761', row[3])

    def test_alpelisib_PIK3CA(self):
        """
        alpelisib	PIK3CA	25544637
        """
        df = self.valid_pki_pk
        item = df.loc[(df['PKI'] == 'alpelisib') & (df['PK'] == 'PIK3CA')]
        self.assertEqual(1, len(item))
        row = item.to_numpy()[0]
        self.assertEqual('n/a', row[2])  # We do not have Kd data for this combination
        self.assertEqual('25544637', row[3])

    def test_apatinib_KIT_removed(self):
        """
        This entry should have been removed since the Kd is too high
        apatinib	KIT	20923544	0.428999964569677	IC50
        """
        df = self.valid_pki_pk
        items = df.loc[(df['PKI'] == 'apatinib') & (df['PK'] == 'KIT')]
        self.assertEqual(0, len(items.index))
