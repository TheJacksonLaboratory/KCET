import os
import pandas as pd
from pathlib import Path
from collections import defaultdict
import datetime
import copy
from typing import List, Dict
from .clinical_trial import ClinicalTrial
from .kinase_inhibitor import KinaseInhibitor


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
        if isinstance(other, Entry):
            return ((self._pki == other._pki) and (self._pk == other._pk) and (self._pmid == other._pmid))
        else:
            return False

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
    A parser for the clinical_trials_data file that is produced by yactp.  (https://github.com/monarch-initiative/yactp).
    """
    def __init__(self,
                 clinical_trials: str,
                 year: int = None):
        """
        :param clinical_trials_data_path: This is the output file from yactp.  (https://github.com/monarch-initiative/yactp).
        :param year: The target year for analysis -- take all data up to this year
        The output file contains kinase-cancer links corresponding to the specific phase that the drugs were tested has
        before that year.
        """
        dir_path = os.path.dirname(os.path.realpath(__file__))  # directory of current file
        basepath = Path(dir_path).parent  # parent directory -- the base of the project
        self.prot_kinase_path = os.path.join(basepath, 'input',
                                             'prot_kinase.tsv')  # file consisting of protein kinases their gene symbols, ncbi gene ids and ensembl gene ids
        self.gene_symbol_to_ncbigene_map = self._parse_prot_kinase()  # map gene symbols to their ncbi gene ids
        self._clinical_trials_data_path = clinical_trials
        if not os.path.isfile(self._clinical_trials_data_path):
            raise FileNotFoundError("Could not find %s" % self._clinical_trials_data_path)
        self.drug_kinase_links_data_path = os.path.join(basepath, 'input', 'drug_kinase_links.tsv')
        if year is None:
            now = datetime.datetime.now()  # default to current year
            self._year = now.year
        else:
            self._year = year
        self._genesymbol_to_id_map = None
        self._drug_kinase_links = None
        self.all_phases_df = None
        self._validation_pki_dict = None  # List of PKI<->disease AFTER target year (can be empty if target year is now)
        self._pki_dict = None  # List of PKI<->disease UP TO target year
        ## Ingest data
        self._ingest_kinase_cancer_links()
        targeted_set, targeted_set_phase_4 = self._ingest_targeted_kinase_set(clinical_trials)
        self._targeted = sorted(list(targeted_set))
        self._targeted_phase_4 = sorted(list(targeted_set_phase_4))

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

    def _parse_kinase_inhibitor_to_kinase_file(self) -> Dict:
        if not os.path.exists(self.drug_kinase_links_data_path):
            raise FileNotFoundError("Cound not find drug_kinase_links.tsv file")
        pki_to_kinase_dict = defaultdict(list)
        with open(self.drug_kinase_links_data_path) as f:
            header = next(f)
            if not header.startswith('PKI'):
                raise ValueError("Malformed header of drug_kinase_links.tsv file")
            for line in f:
                e = Entry(line)
                pki_to_kinase_dict[e.pki].append(e.pk)
        return pki_to_kinase_dict

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
        self._genesymbol_to_id_map = self._parse_prot_kinase()
        self._drug_kinase_links = self._parse_kinase_inhibitor_to_kinase_file()
        trials = self._get_ct_by_phase()
        # sanity check-make a set of the medications
        medications = {t.drug for t in trials}
        print("[INFO] Parsed data for %d medications." % len(medications))
        pki_dict = defaultdict(KinaseInhibitor)
        validation_pki_dict = defaultdict(KinaseInhibitor)
        for t in trials:
            med = t.drug
            if t.start_date <= self._year:
                if not med in pki_dict:
                    ki = KinaseInhibitor(name=med)
                    pki_dict[med] = ki
                ki = pki_dict.get(med)
                ki.add_study(cancer=t.disease, mesh_id=t.mesh_id, nct=t.nct_id, year=t.start_date, phase=t.phase)
            else:
                if not med in validation_pki_dict:
                    ki = KinaseInhibitor(name=med)
                    validation_pki_dict[med] = ki
                ki = validation_pki_dict.get(med)
                ki.add_study(cancer=t.disease, mesh_id=t.mesh_id, nct=t.nct_id, year=t.start_date, phase=t.phase)
        self._pki_dict = pki_dict
        self._validation_pki_dict = validation_pki_dict

    def _get_data_frame(self, dict_list: List, removeRedundantEntries: bool = False):
        """
        Return a pandas dataframe with data for all trials and all phases
        Constructs the dataframe from a list of dictionaries
        """
        extended_dict_list = []  ## The dictionaries where we map the PKIs to kinases/genes ids
        for dct in dict_list:
            medication = dct['pki']
            #print(medication)
            if medication is None:
                raise ValueError("Could not extract PKI")  # should never happen
            if not medication in self._drug_kinase_links:
                # should never happen
                raise ValueError("Could not find " + medication + " in pki to pk dict")
            lst = self._drug_kinase_links[medication]
            for kinase in lst:
                if not kinase in self._genesymbol_to_id_map:
                    # should never happen
                    print(kinase)
                    raise ValueError("Could not find " + kinase + " in gene id map")
                geneid = self._genesymbol_to_id_map[kinase]
                extended_d = copy.deepcopy(dct)
                extended_d['kinase'] = kinase
                extended_d['gene_id'] = geneid
                extended_dict_list.append(extended_d)
        if removeRedundantEntries:
            unique_list = []
            seen_entries = set()
            for entry in extended_dict_list:
                # form a key that is unique for the gene/kinase combination
                key = entry['mesh_id'] + entry['gene_id']
                if not key in seen_entries:
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
        df = self._get_data_frame(dict_list=dict_list, removeRedundantEntries=removeRedundantEntries)
        return df


    def get_phase_4(self, removeRedundantEntries: bool = False):
        """
        Return a pandas dataframe with data for all trials in phase 4
        """
        dict_list = []  ## The dictionaries WITHOUT the info about kinases/gene ids
        for _, v in self._pki_dict.items():
            dict_list.extend(v.get_data_frame_phase_4())
        return self._get_data_frame(dict_list=dict_list, removeRedundantEntries=removeRedundantEntries)

    def get_year(self):
        return self._year

    def get_all_phases_for_training(self):
        df = self.get_all_phases()
        df_tr = df[['gene_id', 'mesh_id']].copy()
        df_tr.drop_duplicates(inplace=True)
        return df_tr

    def get_phase_4_for_training(self):
        df = self.get_phase_4()
        df_tr = df[['gene_id', 'mesh_id']].copy()
        df_tr.drop_duplicates(inplace=True)
        return df_tr

    def get_validation_all_phases(self):
        """
        Get all clinical studies from years after the target year.
        If the target year is now -- throw an error
        """
        if self._year == datetime.datetime.now().year:
            raise ValueError(
                "Cannot get studies from the future! This function can only be used for historical comparisons")
        dict_list = []  ## The dictionaries WITHOUT the info about kinases/gene ids
        for _, v in self._validation_pki_dict.items():
            dict_list.extend(v.get_data_frame_all_phases())
        return self._get_data_frame(dict_list=dict_list)

    def get_validation_phase_4(self):
        """
        Get all clinical studies from years after the target year.
        If the target year is now -- throw an error
        """
        if self._year == datetime.datetime.now().year:
            raise ValueError(
                "Cannot get studies from the future! This function can only be used for historical comparisons")
        dict_list = []  ## The dictionaries WITHOUT the info about kinases/gene ids
        for _, v in self._validation_pki_dict.items():
            dict_list.extend(v.get_data_frame_phase_4())
        return self._get_data_frame(dict_list=dict_list)

    def _ingest_targeted_kinase_set(self, path):
        """
        path -- path to the clinical_trials_by_phase file from yactp
        This function creates two maps -- one with the earliest date by which 
        a 
        """
        # df_all is a Pandas dataframe with all clinical trials
        df_all = self.get_all_phases()
        # not_later_than_target_year = df_all['year']<=self._year
        # df_all_not_later = df_all[not_later_than_target_year]
        targeted_set = set(df_all['kinase'])
        ## Now restrict the targets to phase 4
        is_phase_4 = df_all['phase'] == 'Phase 4'
        df_all_not_later_phase_4 = df_all[is_phase_4]
        targeted_set_phase_4 = set(df_all_not_later_phase_4['kinase'])
        return targeted_set, targeted_set_phase_4

    def get_targeted_kinases_with_gene_id(self):
        dict_list = []
        for kinase in self._targeted:
            gene_id = self._genesymbol_to_id_map.get(kinase, "n/a")
            d = {'kinase': kinase, "gene_id": gene_id}
            dict_list.append(d)
        return pd.DataFrame(dict_list)

    def get_untargeted_kinases_with_gene_id(self):
        dict_list = []
        targeted_set = set(self._targeted)
        for kinase, gene_id in self._genesymbol_to_id_map.items():
            if kinase in targeted_set:
                continue
            d = {'kinase': kinase, "gene_id": gene_id}
            dict_list.append(d)
        return pd.DataFrame(dict_list)

    def get_targeted_kinases_with_gene_id_phase_4(self):
        dict_list = []
        for kinase in self._targeted_phase_4:
            gene_id = self._genesymbol_to_id_map.get(kinase, "n/a")
            d = {'kinase': kinase, "gene_id": gene_id}
            dict_list.append(d)
        return pd.DataFrame(dict_list)
