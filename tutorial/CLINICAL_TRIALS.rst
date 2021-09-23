####################
Clinical Trials data
####################


The purpose of this step is to generate a table of data about clinical trials that investigated the use of protein-kinase inhibitors (PKIs) to treat cancer. For this, we developed a separate Python tool called
`yactp (yet another clinical trials parser) <https://github.com/monarch-initiative/yactp>`_ tool.

Setup the yatcp library as described in the README file of the repository.

The purpose of this script is to retrieve all data in
the `ClinicaTrials.gov <https://clinicaltrials.gov/>`_ resource that describe studies on cancer that are related to our
list of protein-kinase inhibitors (PKIs). We use our list of PKIs that is available
`here <https://github.com/TheJacksonLaboratory/KCET/blob/main/input/protein_kinase_inhibitors.txt>`_.
This file has names of 75 PKIs, one per line: ::

    abemaciclib
    acalabrutinib
    (...)
    vemurafenib
    zanubrutinib


Our goal is to search for clinical trials that describe the use of any of these PKIs to treat a cancer.
ClinicalTrials.gov uses Medical Subject Heading (MeSH) terms to specify diseases.
MeSH is the National Library of Medicine's controlled vocabulary thesaurus.
We first retrieve a list of cancer-related MeSH terms by using the following python script to
generate a file called ``protein_kinase_inhibitors.txt``. We use the ``meshSparql.py`` script
that is available in the yactp repository and run it with the MeSH id for ``Neoplasms``,
`D009369 <https://meshb.nlm.nih.gov/record/ui?ui=D009369>`_. The script returns all of the
specific descendents of this term.  ::

    python meshSparql.py -m D009369



This retrieves 697 MeSH entries that descend from Neoplasm, writing them to a file with line such
as the following. ::

    000182	ACTH Syndrome, Ectopic
    D049913	ACTH-Secreting Pituitary Adenoma
    D000008	Abdominal Neoplasms
    D058739	Aberrant Crypt Foci
    D049309	Acanthoma
    D018250	Acrospiroma
    D050398	Adamantinoma
    D000230	Adenocarcinoma



The file is called  \verb+mesh_D009369_MM_DD_YYYY.tsv+ (MM, DD, and YYYY will be replaced by the date on which the script was run).

Finally, we run the following script. ::

    python parseCT.py -i protein_kinase_inhibitors.txt -m mesh_D009369_MM_DD_YYYY.tsv

This generates the file ``clinical_trials_by_phase.tsv`` that is used in subsequent steps.
The version of this file that we used for the analysis in the manuscript is available in the
project zenodo repository (see README).