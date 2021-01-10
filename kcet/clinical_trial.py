

class ClinicalTrial:
    """
    Represents one line from the clinical_trials_by_phase.tsv file
    disease	mesh_id	drug	phase	start_date	completion_date	nct_id
    Carcinoma, Non-Small-Cell Lung	D002289	afatinib	Phase 1	2009	2020	NCT01999985;NCT03054038;NCT04448379;NCT02191891;NCT03827070;NCT01090011;NCT03711422;NCT00993499;NCT01647711;NCT02364609;NCT01288430
    """
    def __init__(self, 
                disease:str,	
                mesh_id:str,
                drug:str,	
                phase:str,
                start_date:int,	
                completion_date:int,
                nct_id:str):
        self._disease = disease
        self._mesh_id = mesh_id
        self._drug = drug
        allowed_values = {'Phase 1', 'Phase 2', 'Phase 3', 'Phase 4'}
        if not phase in allowed_values:
            raise ValueError("Bad value for phase -- '%s'" % phase)
        self._phase = phase
        self._start_date = start_date
        self._completion_data = completion_date
        self._nct_id = nct_id

    @property
    def disease(self):
        return self._disease

    @property
    def mesh_id(self):
        return self._mesh_id

    @property
    def drug(self):
        return self._drug

    @property
    def phase(self):
        return self._phase

    @property
    def start_date(self):
        return self._start_date

    @property
    def completion_date(self):
        return self._completion_data

    @property
    def nct_id(self):
        return self._nct_id

    def is_phase_4(self):
        return self._phase == 'Phase 4'

    
