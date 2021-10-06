from .clinical_trial import ClinicalTrial
from .kinase_inhibitor import KinaseInhibitor
from .drugcentral_pk_pki_parser import DrugCentralPkPkiParser

import os
import pandas as pd
from pathlib import Path
from collections import defaultdict
import copy
from typing import List, Dict
import logging

logging.basicConfig(format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
                    datefmt='%Y-%m-%d:%H:%M:%S',
                    filename='kcet.log',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)


class Entry:
    """
    Simple class representing one line in the drug_kinase_links.tsv file
    Attributes:
        _pki     A protein-kinase inhibitor (PKI).
        _pk  A protein kinase that is inhibited by the PKI.
        _act_value  activity (less than 30 considered active)
        _pmid   A PubMed id for the PKI-PK link
    """

    def __init__(self, line: str) -> None:
        super().__init__()
        fields = line.rstrip().split('\t')
        if len(fields) != 4:
            raise ValueError(line)
        self._pki = fields[0]
        self._pk = fields[1]
        self._act_value = fields[2]
        self._pmid = fields[3]

    def __eq__(self, other):
        return isinstance(other,
                          Entry) and self._pki == other._pki and self._pk == other._pk and self._pmid == other._pmid

    def __hash__(self):
        return hash((self._pki, self._pk, self._act_value, self._pmid))

    @property
    def pki(self) -> str:
        return self._pki

    @property
    def pk(self) -> str:
        return self._pk


class CTParserByPhase:
    """
    A parser for the clinical_trials_data file that is produced by yactp.
    (https://github.com/monarch-initiative/yactp).

    """

    def __init__(self, clinical_trials: str, n_pk: int = 5):
        """
        :param clinical_trials: This is the output file from yactp.  (https://github.com/monarch-initiative/yactp).
        The output file contains kinase-cancer links corresponding to the specific phase that the drugs were tested
        """
        dir_path = os.path.dirname(os.path.realpath(__file__))  # directory of current file
        base_path = Path(dir_path).parent  # parent directory -- the base of the project
        # file consisting of protein kinases their gene symbols, ncbi gene ids and ensembl gene ids
        self.prot_kinase_path = os.path.join(base_path, 'input', 'prot_kinase.tsv')
        self.gene_symbol_to_ncbigene_map = self._parse_prot_kinase()  # map gene symbols to their ncbi gene ids
        self._clinical_trials_data_path = clinical_trials
        if not os.path.isfile(self._clinical_trials_data_path):
            raise FileNotFoundError("Could not find %s" % self._clinical_trials_data_path)
        self.drug_kinase_links_data_path = os.path.join(base_path, 'input', 'drug_kinase_links.tsv')
        # _genesymbol_to_id map has entries like {CDK20:ncbigene23552}
        self._genesymbol_to_id_map = self._parse_prot_kinase()

        self.all_phases_df = None
        drug_central_parser = DrugCentralPkPkiParser()
        # both of the following variables have entries like {abemaciclib:[CDK4,CDK6]}
        # the thresholded version is used for training.
        # the 'all links' version has all links from DrugCentral, regardless of threshold -- we do not want
        # to test against this because it would not be fair.
        self._thresholded_pk_pki_links = drug_central_parser.get_pk_pki_with_threshold(n_pki_limit=n_pk)
        self._all_pk_pki_links = drug_central_parser.get_all_pk_pki()
        self._ingest_kinase_cancer_links()

    def _parse_prot_kinase(self) -> Dict:
        """
        parse prot_kinase.tsv and map gene symbols to ncbigene ids
        return a dictionary with key being a gene symbol and value being the corresponding NCBI gene id
        e.g., for the line
        CDK20	cyclin dependent kinase 20	23552	ENSG00000156345
        we would have {'CDK20' => 23552}
        """
        symbol2ncbigene = defaultdict(str)
        with open(self.prot_kinase_path) as f:
            # Note this file has no header
            for line in f:
                fields = line.rstrip().split('\t')
                if len(fields) != 4:
                    # should never happen
                    raise ValueError("Malformed line: %s" % line)
                sym = fields[0].strip()
                ncbigene = fields[2].strip()
                ncbi_gene_id = "ncbigene" + str(
                    ncbigene)  # We add "ncbigene" because in pubmed_cr.tsv, genes are represented in this way
                symbol2ncbigene[sym] = ncbi_gene_id
        return symbol2ncbigene

    def _get_ct_by_phase(self) -> List:
        """
        parse the clinical_trials_by_phase.tsv file
        """
        if not os.path.exists(self._clinical_trials_data_path):
            raise FileNotFoundError("Could not find the clinical_trials_by_phase.tsv file")
        trials = []
        with open(self._clinical_trials_data_path) as f:
            header = next(f)
            # should be disease	mesh_id	drug	phase	start_date	completion_date	nct_id
            if not header.startswith("disease"):
                raise ValueError("Bad header line in clinical_trials_by_phase.tsv file")
            for line in f:
                fields = line.rstrip().split('\t')
                if len(fields) != 7:
                    raise ValueError("Bad line in clinical_trials_by_phase.tsv file:" + line)
                disease = fields[0]
                mesh = fields[1]
                mesh_id_first_letter = mesh[0].lower()
                mesh_id = "mesh" + mesh_id_first_letter + mesh[1:]
                drug = fields[2]
                phase = fields[3]
                start = int(fields[4])
                end = int(fields[5])
                nct = fields[6]
                ct = ClinicalTrial(disease=disease, mesh_id=mesh_id, drug=drug, phase=phase,
                                   start_date=start, completion_date=end, nct_id=nct)
                trials.append(ct)
        return trials

    def _ingest_kinase_cancer_links(self) -> None:
        """
        We first use the clinical_trials_by_phase.tsv to extract disease-drug information for a specific phase,
        e.g. Phase 1 before a specific date (year). The phase and date are specified by the user.
        Then, we use the drug_kinase_links.tsv to obtain drug-kinase links and finally we generate disease_kinase links.
        """
        trials = self._get_ct_by_phase()
        # sanity check-make a set of the medications
        medications = {t.drug for t in trials}
        logging.info("Parsed data for %d medications." % len(medications))
        pki_dict = defaultdict(KinaseInhibitor)
        for t in trials:
            med = t.drug
            if med not in pki_dict:
                ki = KinaseInhibitor(name=med)
                pki_dict[med] = ki
            ki = pki_dict.get(med)
            ki.add_study(cancer=t.disease, mesh_id=t.mesh_id, nct=t.nct_id, year=t.start_date, phase=t.phase)
        self._pki_dict = pki_dict

    def _get_data_frame(self, dict_list: List, remove_redundant_entries: bool = False) -> pd.DataFrame:
        """
        Return a pandas dataframe with data for all trials and all phases
        Constructs the dataframe from a list of dictionaries
        """
        extended_dict_list = []  ## The dictionaries where we map the PKIs to kinases/genes ids
        for dct in dict_list:
            medication = dct['pki']
            if medication is None:
                raise ValueError("Could not extract PKI")  # should never happen
            if medication not in self._thresholded_pk_pki_links:
                # should never happen
                raise ValueError("Could not find " + medication + " in pki to pk dict")
            lst = self._thresholded_pk_pki_links[medication]
            for kinase in lst:
                if kinase not in self._genesymbol_to_id_map:
                    # should never happen
                    raise ValueError("Could not find " + kinase + " in gene id map")
                gene_id = self._genesymbol_to_id_map[kinase]
                extended_d = copy.deepcopy(dct)
                extended_d['kinase'] = kinase
                extended_d['gene_id'] = gene_id
                extended_dict_list.append(extended_d)
        if remove_redundant_entries:
            unique_list = []
            seen_entries = set()
            for entry in extended_dict_list:
                # form a key that is unique for the gene/kinase combination
                key = entry['mesh_id'] + entry['gene_id']
                if key not in seen_entries:
                    seen_entries.add(key)
                    unique_list.append(entry)
            df = pd.DataFrame.from_records([d for d in unique_list])
        else:
            df = pd.DataFrame.from_records([d for d in extended_dict_list])
        # reorder the columns
        newcols = ['cancer', 'mesh_id', 'kinase', 'gene_id', 'pki', 'nct', 'phase', 'year']
        return df[newcols]

    def get_all_phases(self, removeRedundantEntries: bool = False):
        """
        Return a pandas dataframe with data for all trials and all phases
        """
        dict_list = []  ## The dictionaries WITHOUT the info about kinases/gene ids
        for _, v in self._pki_dict.items():
            dict_list.extend(v.get_data_frame_all_phases())
        df = self._get_data_frame(dict_list=dict_list, remove_redundant_entries=removeRedundantEntries)
        return df

    def get_phase_4(self, removeRedundantEntries: bool = False):
        """
        Return a pandas dataframe with data for all trials in phase 4
        """
        dict_list = []  ## The dictionaries WITHOUT the info about kinases/gene ids
        for _, v in self._pki_dict.items():
            dict_list.extend(v.get_data_frame_phase_4())
        return self._get_data_frame(dict_list=dict_list, remove_redundant_entries=removeRedundantEntries)

    def get_year(self):
        return self._year
