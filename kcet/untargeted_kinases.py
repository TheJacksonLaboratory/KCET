import pandas as pd
from collections import defaultdict
import os
import datetime


from .ct_by_phase_parser import CTParserByPhase
from .protein_kinase_parser import ProteinKinaseParser

class UnTargetedKinases:
    def __init__(self, clinical_trials_by_phase: str, year: int = None):
        """
        Extract protein kinases that have not been targeted.
        The assumption of this script is that all targeted kinases
        are listed in the file drug_kinase_links.tsv:
        abemaciclib	CDK4	24919854
        abemaciclib	CDK6	24919854
        The kinases that are listed in prot_kinase.tsv but NOT in drug_kinase_links.tsv are
        currently UNTARGETED.
        The script is hard-wired to input the files in the ``input`` subdirectory.
        """
        if not year is None:
            self._year = int(year)
        else:
            now = datetime.datetime.now()
            year = int(now.year)
            self._year = year
        # This file is in the kcet directory. We need to get to the parent directory
        # The following gets a list of all protein kinases
        pkparser = ProteinKinaseParser()
        
        self._genesymbol_to_id_map = pkparser.get_symbol_to_id_map()
        targeted_set, targeted_set_phase_4 = self._get_targeted_kinase_set(clinical_trials_by_phase)
        self._targeted = sorted(list(targeted_set))
        self._targeted_phase_4 = sorted(list(targeted_set_phase_4))


    

    def _get_targeted_kinase_set(self, path):
        """
        path -- path to the clinical_trials_by_phase file from yactp
        This function creates two maps -- one with the earliest date by which 
        a 
        """
        parser = CTParserByPhase(clinical_trials=path)
        # df_all is a Pandas dataframe with all clinical trials
        df_all = parser.get_all_phases()
        not_later_than_target_year = df_all['year']<=self._year
        df_all_not_later = df_all[not_later_than_target_year]
        targeted_set = set(df_all_not_later['kinase'])
        ## Now restrict the targets to phase 4
        is_phase_4 = df_all_not_later['phase'] == 'Phase 4'
        df_all_not_later_phase_4 = df_all_not_later[is_phase_4]
        targeted_set_phase_4 = set(df_all_not_later_phase_4['kinase'])
        return targeted_set, targeted_set_phase_4

    def get_targeted_kinases_with_gene_id(self):
        dict_list = []
        for kinase in self._targeted:
            gene_id = self._genesymbol_to_id_map.get(kinase, "n/a")
            d = { 'kinase': kinase, "gene.id": gene_id}
            dict_list.append(d)
        return pd.DataFrame(dict_list)

    def get_untargeted_kinases_with_gene_id(self):
        dict_list = []
        targeted_set = set(self._targeted)
        for kinase, gene_id in self._genesymbol_to_id_map.items():
            if kinase in targeted_set:
                continue
            d = { 'kinase': kinase, "gene.id": gene_id}
            dict_list.append(d)
        return pd.DataFrame(dict_list)

    def get_targeted_kinases_with_gene_id_phase_4(self):
        dict_list = []
        for kinase in self._targeted_phase_4:
            gene_id = self._genesymbol_to_id_map.get(kinase, "n/a")
            d = { 'kinase': kinase, "gene.id": gene_id}
            dict_list.append(d)
        return pd.DataFrame(dict_list)

    def get_target_year(self):
        return self._year


