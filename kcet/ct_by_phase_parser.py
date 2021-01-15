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
    """
    def __init__(self, line: str) -> None:
        super().__init__()
        fields = line.rstrip().split('\t')
        if len(fields) != 3:
            raise ValueError(line)
        self._pki = fields[0]
        self._pk = fields[1]
        self._pmid = fields[2]

    def __eq__(self, other):
        if isinstance(other, Entry):
            return ((self._pki == other._pki) and (self._pk == other._pk) and (self._pmid == other._pmid))
        else:
            return False

    def __hash__(self):
        return hash((self._pki, self._pk, self._pmid))

    @property
    def pki(self) -> str:
        return self._pki

    @property
    def pk(self) -> str:
        return self._pk

class CTParserByPhase:
    def __init__(self, 
                clinical_trials: str, 
                year: int = None):
        """
        :param clinical_trials_data_path: This is the output file from yactp.  (https://github.com/monarch-initiative/yactp).
        :param drug_kinase_links_data_path: This is a curated list of drug and kinase links.
        :param date: The date( year) that the drug was tested on diseases.
        The output file contains kinase-cancer links corresponding to the specific phase that the drugs were tested has
        before that year.
        """
        dir_path = os.path.dirname(os.path.realpath(__file__))  # directory of current file
        basepath = Path(dir_path).parent  # parent directory -- the base of the project
        self.prot_kinase_path = os.path.join(basepath, 'input', 'prot_kinase.tsv') #file consisting of protein kinases their gene symbols, ncbi gene ids and ensembl gene ids
        self.gene_symbol_to_ncbigene_map = self._parse_prot_kinase() #map gene symbols to their ncbi gene ids

        self.clinical_trials_data_path = clinical_trials
        self.drug_kinase_links_data_path = os.path.join(basepath, 'input', 'drug_kinase_links.tsv')
        if year is None:
            now = datetime.datetime.now() # default to current year
            self._year = now.year
        else:
            self._year = year
        self._genesymbol_to_id_map = None
        self._drug_kinase_links = None
        self.all_phases_df = None
        ## Ingest data
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
                ncbi_gene_id = "ncbigene" + str(ncbigene) #We add "ncbigene" because in pubmed_cr.tsv, genes are represented in this way
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
        if not os.path.exists(self.clinical_trials_data_path):
            raise FileNotFoundError("Could not find the clinical_trials_by_phase.tsv file")
        trials = []
        with open(self.clinical_trials_data_path) as f:
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
                ct = ClinicalTrial(disease=disease, mesh_id=mesh_id, drug=drug,phase=phase,
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
        print("[INFO] Parsed data for the following medications:")
        for med in medications:
            print("[INFO] %s" % med)
        pki_dict = defaultdict(KinaseInhibitor)
        for t in trials:
            if t.start_date > self._year:
                continue  # Skip study newer than target year.
            med = t.drug
            if not med in pki_dict:
                ki = KinaseInhibitor(name=med)
                pki_dict[med] = ki
            ki = pki_dict.get(med)
            ki.add_study(cancer=t.disease, mesh_id=t.mesh_id, nct=t.nct_id, year=t.start_date, phase=t.phase)
        self._pki_dict = pki_dict

    def _get_data_frame(self, dict_list: List):
        """
        Return a pandas dataframe with data for all trials and all phases
        Constructs the dataframe from a list of dictionaries
        """   
        extended_dict_list = [] ## The dictionaries where we map the PKIs to kinases/genes ids
        for dct in dict_list:
            medication = dct['pki']
            if medication is None:
                raise ValueError("Could not extract PKI") # should never happen
            if not medication in self._drug_kinase_links:
                # should never happen
                raise ValueError("Could not find " + medication + " in pki to pk dict")
            lst = self._drug_kinase_links[medication]
            for kinase in lst:
                if not kinase in self._genesymbol_to_id_map:
                    # should never happen
                    raise ValueError("Could not find " + kinase + " in gene id map")
                geneid = self._genesymbol_to_id_map[kinase]
                extended_d = copy.deepcopy(dct)
                extended_d['kinase'] = kinase
                extended_d['gene.id'] = geneid
                extended_dict_list.append(extended_d)
        df = pd.DataFrame.from_records([d for d in extended_dict_list])
        # reorder the columns
        newcols = ['cancer', 'mesh_id', 'kinase', 'gene.id',  'pki', 'nct', 'phase', 'year']
        return df[newcols]
    
 
    def get_all_phases(self):
        """
        Return a pandas dataframe with data for all trials and all phases
        """  
        dict_list = [] ## The dictionaries WITHOUT the info about kinases/gene ids
        for _, v in self._pki_dict.items():
            dict_list.extend(v.get_data_frame_all_phases()) 
        return self._get_data_frame(dict_list=dict_list)

    def get_phase_4(self):
        """
        Return a pandas dataframe with data for all trials in phase 4
        """   
        dict_list = [] ## The dictionaries WITHOUT the info about kinases/gene ids
        for _, v in self._pki_dict.items():
            dict_list.extend(v.get_data_frame_phase_4())
        return self._get_data_frame(dict_list=dict_list)

    def get_year(self):
        return self._year

    def get_all_phases_for_training(self):
        df = self.get_all_phases()
        df_tr = df[['gene.id', 'mesh_id']].copy()
        df_tr.drop_duplicates(inplace=True)
        return df_tr

    def get_phase_4_for_training(self):
        df = self.get_phase_4()
        df_tr = df[['gene.id', 'mesh_id']].copy()
        df_tr.drop_duplicates(inplace=True)
        return df_tr

