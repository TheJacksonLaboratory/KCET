from collections import defaultdict
import os


class NeoplasmParser:

    def __init__(self) -> None:
        d = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) 
        self._neoplasms_labels_tsv_path = os.path.join(d, 'input', 'neoplasms_labels.tsv')
        if not os.path.exists(self._neoplasms_labels_tsv_path):
            raise FileNotFoundError("Could not find file at %s" % self._neoplasms_labels_tsv_path)
        print("[INFO] Reading protein kinase information from %s" % self._neoplasms_labels_tsv_path)


    def get_mesh_id_list(self):
        """
        Return a list of the MeSH ids that correspond to cancers from the 
        inputs/neoplasms_labels,tsv file
        """
        mesh_list = []
        with open(self._neoplasms_labels_tsv_path) as f:
            for line in f:
                fields = line.rstrip().split('\t')
                if len(fields) < 2 or len(fields)> 3:
                    #At a minimum, there is MeSH id and label. Most lines have a third field (synonyms)
                    raise ValueError("Bad line in neoplasms_labels,tsv file: %s"  % line)
                mesh = fields[0]
                mesh_id_first_letter = mesh[0].lower()
                mesh_id = "mesh" + mesh_id_first_letter + mesh[1:]
                mesh_list.append(mesh_id)
        return mesh_list

