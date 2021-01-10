from collections import defaultdict


class Cancer:
    def __init__(self, name:str, mesh_id:str) -> None:
        self._name = name
        self._mesh_id = mesh_id
        self._studies_by_phase = defaultdict(list)

    def add_study(self, nct:str, phase:str, year:int) -> None:
        s = Study(nct=nct, phase=phase,year=year)
        self._studies_by_phase[phase].append(s)

    @property
    def name(self):
        return self._name

    @property
    def mesh_id(self):
        return self._mesh_id

    def get_studies_dict(self):
        return self._studies_by_phase


class Study:
    def __init__(self, nct, phase, year) -> None:
        
        self._nct = nct
        self._phase = phase
        self._year = year # start year

    @property
    def nct(self):
        return self._nct

    @property
    def phase(self):
        return self._phase

    @property
    def year(self):
        return self._year

    def is_phase_4(self):
        return self._phase == "Phase 4"


class KinaseInhibitor:
    """
    Represents one kinase inhibitor as well as the cancers and phases of all associated clinical trials.
    """
    def __init__(self, name:str) -> None:
        self._name = name
        self._cancerdict = defaultdict(Cancer)

    def add_study(self, cancer:str, mesh_id:str, nct:str, phase:str, year:int) -> None:
        if not mesh_id in self._cancerdict:
            c = Cancer(name=cancer, mesh_id=mesh_id)
            self._cancerdict[mesh_id] = c
        c = self._cancerdict.get(mesh_id)
        c.add_study(nct=nct, phase=phase,year=year)
        

    def get_data_frame_all_phases(self):
        list_of_dicts = []
        pki = self._name
        for _, v in self._cancerdict.items():
            cancer = v.name
            mesh = v.mesh_id
            study_d = v.get_studies_dict()
            for _, study_list in study_d.items():
                for study in study_list:
                    d = {
                    'pki': pki,
                    'cancer': cancer,
                    'mesh_id': mesh,
                    'phase': study.phase,
                    'year': study.year,
                    'nct': study.nct
                    }
                    list_of_dicts.append(d)
        return list_of_dicts

    def get_data_frame_phase_4(self):
        list_of_dicts = []
        pki = self._name
        for _, v in self._cancerdict.items():
            cancer = v.name
            mesh = v.mesh_id
            study_d = v.get_studies_dict()
            for _, study_list in study_d.items():
                for study in study_list:
                    if study.is_phase_4():
                        d = {
                            'pki': pki,
                            'cancer': cancer,
                            'mesh_id': mesh,
                            'phase': study.phase,
                            'year': study.year,
                            'nct': study.nct
                        }
                        list_of_dicts.append(d)
        return list_of_dicts



        