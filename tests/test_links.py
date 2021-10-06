from kcet.kcet_dataset_generator import Link
import os
from unittest import TestCase


class TestLinks(TestCase):
    """
    Confirm that equality/hashing are working as expected
    """

    @classmethod
    def setUpClass(cls):
        current_dir = os.path.dirname(__file__)
        ct_by_phase_path = os.path.join(current_dir, 'data', "small_ct_by_phase.tsv")
        target_year = 2020
        link1 = Link(cancer='meshd016411', kinase='ncbigene9263')
        link2 = Link(cancer='meshd016066', kinase='ncbigene780')
        link3 = Link(cancer='meshd005870', kinase='ncbigene3932')
        link4 = Link(cancer='meshd012208', kinase='ncbigene27')
        link5 = Link(cancer='meshd002295', kinase='ncbigene8767')
        cls._link_set = set()
        cls._link_set.add(link1)
        cls._link_set.add(link2)
        cls._link_set.add(link3)
        cls._link_set.add(link4)
        cls._link_set.add(link5)

    def test_five_links(self):
        self.assertEqual(5, len(self._link_set))

    def test_simple_equality(self):
        link1a = Link(cancer='meshd016411', kinase='ncbigene9263')
        link1b = Link(cancer='meshd016411', kinase='ncbigene9263')
        link2 = Link(cancer='meshd016066', kinase='ncbigene780')
        self.assertEqual(link1a, link1a)
        self.assertEqual(link1a, link1b)
        self.assertEqual(link1b, link1b)
        self.assertNotEqual(link1a, link2)
        self.assertNotEqual(link1b, link2)
        self.assertEqual(link2, link2)

    def test_exists_in_set(self):
        link1 = Link(cancer='meshd016411', kinase='ncbigene9263')
        link2 = Link(cancer='meshd016066', kinase='ncbigene780')
        link3 = Link(cancer='meshd0FAKE', kinase='ncbigene9263')
        link4 = Link(cancer='meshd016066', kinase='ncbigene78FAKE')
        self.assertTrue(link1 in self._link_set)
        self.assertTrue(link2 in self._link_set)
        self.assertFalse(link3 in self._link_set)
        self.assertFalse(link4 in self._link_set)





