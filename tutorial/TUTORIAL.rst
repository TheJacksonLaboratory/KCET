########
Tutorial
########

This document intends to explain how the analysis was performed in  Ravanmehr et al (2021) ``Supervised learning with word embeddings derived from PubMed captures latent knowledge about protein kinases and cancer`` and to suggest how our approach can be extended to other areas.
The analysis goes through several steps.


zenodo repository
#################

Steps 1 and 2 in our analysis are compute-intensive or require substantial work or computation to (re)produce. We therefore provide
a `zenodo repository <https://zenodo.org/record/5329035>`_ with the files that are produced by these steps.


1. Concept replacement and PubMed abstract processing
#####################################################

The first step of our analysis consists in the creation of concept embeddings based on ``relevant`` abstracts in PubMed.
PubMed abstracts were downloaded (see `here <https://www.nlm.nih.gov/databases/download/pubmed_medline.html>`_ for instructions)
and processed to replace words and phrases by corresponding concept ids based on the `PubTator resource <https://www.ncbi.nlm.nih.gov/research/pubtator/>`_. 

The methods followed in  Ravanmehr et al are described in the Methods section of the manuscript. 
For the experiments described in the main text, we used an in-house Python package called
\href{https://github.com/TheJacksonLaboratory/marea}{marea} to leverage PubTator resources to perform concept replacement and additionally to perform stop-word removal and to
format the abstracts in preparation for word/concept embedding. Since the submission of this paper, 
PubTator has adopted a new format (BioCXML) and the marea scripts will no longer work. Users should adapt concept replacement to the needs of their analysis goals, and can use the marea scripts as a jumping off point.

2. Concept embedding
####################


Word embedding was performed with our `embiggen package <https://pypi.org/project/embiggen/>`_ , which leverages
tensorflow 2 for a variety of machine learning algorithms. This `page <CONCEPT_EMBEDDING.rst>`_ presents the script together with explanations.

3. Clinical Trials data
#######################

The purpose of this step is to generate a table of data about clinical trials that investigated the use of protein-kinase inhibitors (PKIs) to treat cancer.
This `page <CLINICAL_TRIALS.rst>`_ presents details.

4. Random forest classification
###############################

The final step of the pipeline performs random forest classification on vectors representing protein-kinase cancer pairs.
The notebook `randomForestClassification <../notebooks/randomForestClassification.ipynb>`_ explains this step in detail.
See also the script `runRandomForest <../scripts/runRandomForest.py>`_. Novel predictions (2021 onwards) were generated with
code that is explained in the notebook `novelPredictions <../notebooks/novelPredictions.ipynb>`_.


