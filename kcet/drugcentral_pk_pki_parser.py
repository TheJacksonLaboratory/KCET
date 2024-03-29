import os
import csv
from collections import defaultdict
from typing import List, Dict
import pandas as pd
import logging

logging.basicConfig(format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
                    datefmt='%Y-%m-%d:%H:%M:%S',
                    filename='kcet.log',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)


class PkPkiLink:
    """
    This class models a protein kinase (PKs) to protein kinase inhibitor (PKI) link
    Attributes:
        _pki     A protein-kinase inhibitor (PKI).
        _pk  A protein kinase that is inhibited by the PKI.
        _act_value  activity (less than 30 considered active)
        _pmid   A PubMed id for the PKI-PK link
    """

    def __init__(self, pki: str, pk: str, act_val: float, pmid: str):
        self._pki = pki
        self._pk = pk
        self._act_value = act_val
        self._pmid = pmid

    @property
    def pk(self):
        return self._pk

    @property
    def pki(self):
        return self._pki

    @property
    def act_value(self):
        return self._act_value

    @property
    def pmid(self):
        return self._pmid

    def is_valid(self, act_threshold: float = 0.03):
        """
        Do we have an affinity value and is it less than the threshold
        """
        if self._act_value is None:
            return False
        return self._act_value <= act_threshold

    def is_valid_or_not_provided(self, act_threshold: float = 0.03):
        """
        Note that we assume a PK-PKI can be valid if it does not have any value for
        self._act_value, because this means the link came from PubMed rather than 
        the DrugCentral affinity data. However, other code will favor PK-PKI links
        that have explicit affinity data
        """
        if self._act_value is None:
            return True
        else:
            return self._act_value <= act_threshold

    def act_is_none(self):
        return self._act_value is None

    def __lt__(self, other):
        if other._act_value is None:
            return True
        if self._act_value is None:
            return False
        return self.act_value < other.act_value

    def __eq__(self, other):
        return (self._pki, self._pk, self._act_value, self._pmid) == (
        other._pki, other._pk, other._act_value, other._pmid)

    def __str__(self) -> str:
        return "%s-%s (%s)" % (self._pki, self._pk, str(self._act_value))

    @property
    def tsv_row(self):
        items = [self._pki, self._pk, str(self._act_value), self._pmid]
        return "\t".join(items)


class DrugCentralPkPkiParser:
    """
    The purpose of this class is to determine the list of protein kinases (PKs) that are closely associated to
    a given protein kinase inhibitor (PKI) based on knowledge about their dissociation constant (Kd).
    We use 30 nM as a filter, keeping all PKI<->PK links where the Kd is not higher than 30nM.
    And, for drugs with no PKs below the cutoff, let's relax the cutoff in 1 log10 increments, so 300 nM and if still
    no kinases, go to 3 uM.  I doubt that would be the case - and I would certainly not see the value in the 3 uM scenario.
    We use a file with Kds from DrugCentral.
    Note that the Act_Value is micromolar, and thus, the appropriate cut-off is 0.03, i.e., 0.03 micromolar = 30 nanomolar.
    """

    def __init__(self):
        """
        Note that we assume that the DrugCentral file is in the input directory.
        """
        current_dir = os.path.dirname(__file__)
        parent_dir = os.path.dirname(os.path.normpath(current_dir))
        drug_central_pk_pki_file = os.path.join(parent_dir, 'input', "DrugCentralPKIPK.csv")
        self._pk_pki_list = []
        logger.info("Reading PK/PKI data from %s", drug_central_pk_pki_file)
        with open(drug_central_pk_pki_file) as f:
            rdr = csv.DictReader(f, delimiter='\t')
            for row in rdr:
                pki = row['PKI']
                pk = row['PK']
                pmid = row['PMID']
                act_value = row['ACT_VALUE (uM)']
                if act_value == '':
                    act_value = None
                else:
                    act_value = float(row['ACT_VALUE (uM)'])
                act_type = row['ACT_TYPE']
                if act_type is not None and len(act_type) > 1:
                    if act_type != 'Kd' and act_type != 'Ki' and act_type != 'IC50' and act_type != 'EC50':
                        raise ValueError("Unrecognized ACT_TYPE: %s" % act_type)
                pk_pki = PkPkiLink(pki=pki, pk=pk, act_val=act_value, pmid=pmid)
                self._pk_pki_list.append(pk_pki)
        logger.info("Ingested %d pk pki links with Kd data", len(self._pk_pki_list))

    def _get_max_affinity_links(self, pk_pki, n_pki_limit: int) -> List:
        if not isinstance(pk_pki, list):
            raise ValueError("pk_pki must be a list")
        if len(pk_pki) == 0:
            raise ValueError("pk_pki cannot be empty")
        if len(pk_pki) <= n_pki_limit:
            # no problem, there are not too many PKIs
            return pk_pki
        # If we get here, then there are more PKs than desired at the indicated threshold
        # We want to reduce the Kd threshold to get a smaller number of higher affinity PKI<->PK links
        sorted_pk_pki = sorted(pk_pki)
        return sorted_pk_pki[:n_pki_limit]

    def get_pk_pki_with_threshold(self, n_pki_limit: int = 5, threshold: float = 0.03) -> Dict:
        """
        n_pki_limit: Limit on the number of PKs that are inhibited per PKI
        We do not want to include PKIs that are very promiscuous because the relation between
        a PK and a cancer treated by that PKI becomes tenuous
        Return a map with entries like this {abemaciclib:[CDK4,CDK6]}
        """
        valid_pki_dict = defaultdict(list)
        for pk_pki in self._pk_pki_list:
            if pk_pki.is_valid_or_not_provided(act_threshold=threshold):
                valid_pki_dict[pk_pki.pki].append(pk_pki)
        valid_pk_pki = []
        thresholded_pk_pki = defaultdict(list)
        for k, v in valid_pki_dict.items():
            logger.info("Processing %s", k)
            if len(v) > n_pki_limit:
                logger.warning("Adjusting threshold for PKI {} because it inhibits too many PKs ({})".format(k, len(v)))
                v = self._get_max_affinity_links(pk_pki=v, n_pki_limit=n_pki_limit)
            thresholded_pk_pki[k] = [item.pk for item in v]
            logger.info("PKI: %s, number of inhibited PKs %d", k, len(v))
        logger.info("Extracted %d PKI<->PK interactions", len(valid_pk_pki))
        return thresholded_pk_pki

    def get_all_pk_pki(self, threshold: float = 0.03) -> Dict:
        """
        For the novel predictions, we do not want to filter the PK/PKI links from drug central since we
        are interested in finding things that are not previously investigated, whether the affinity is
        within the n_pki_limit used for training/historical validation or not
        Return a map with entries like this {abemaciclib:[CDK4,CDK6]}
        """
        valid_pki_dict = defaultdict(list)
        for pk_pki in self._pk_pki_list:
            if pk_pki.is_valid_or_not_provided(act_threshold=threshold):
                valid_pki_dict[pk_pki.pki].append(pk_pki.pk)
        return valid_pki_dict


