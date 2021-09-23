# kinaseCancerEmbeddingTool

This repository presents the analysis and code used for 
[Ravanmehr et al. (2021) Supervised learning with word embeddings derived from PubMed captures latent knowledge about protein kinases and cancer](https://www.biorxiv.org/content/10.1101/2021.06.11.447943v1).


### zenodo

The analysis encompasses a series of steps, some of which require
substantial work or computation to (re)produce. Therefore, we have
placed several files that represent intermediate steps of the analysis in
this [zenodo repository](https://zenodo.org/record/5516252).

- clinical_trials_by_phase.tsv
- embedding_SG_dim100_upto2010.npy
- embedding_SG_dim100_upto2014.npy	  
- embedding_SG_dim100_upto2020.npy	  
- words_SG_upto2010.txt
- words_SG_upto2014.txt	  
- words_SG_upto2020.txt

The clinical_trials_by_phase.tsv file represents the output of the YATCP tool (See below).
The embedding and words files represent word/concept embedded vectors and labels (See
the [tensorflow document](https://www.tensorflow.org/text/guide/word_embeddings) for details about 
file formats).

## Start to finish tutorial

See the [tutorial](TUTORIAL.rst) for a streamlined, start to finish explanation of the entire pipeline.

The notebooks in this repository demonstrate the analysis steps that follow generation of word/concept embedddings.


## Jupyter notebooks

Set up your jupyter environment with a kernel with all required packages.
There are many ways of doing this. The following is one way.

```
virtualenv mykernel
source mykernel/bin/activate
(mykernel) $ pip install -r requirments
(mykernel) $ pip install jupyter
(mykernel) $ ipython kernel install --name "local-venv" --user
(mykernel) $ jupyter-lab
```
This will create a kernel called ``local-venv`` that will be visible in the
jupyter lab environment (any name can be used).

We provide the following notebooks:

- [figure1](notebooks/figure1.ipynb): This notebook explores and visualizes the input data and shows how Figure 1C and 1D were generated.
- [randomForestClassification](notebooks/randomForestClassification.ipynb): This notebook shows how to perform random forest classification using our analysis pipeline, and demonstrates how the AUC and PR plots in the manuscript and supplement were generated.
- [novelPredictions](notebooks/novelPredictions.ipynb): This notebook shows how we used a random forest model to generate de novo predictions.

## Scripts
The scripts can be run with the same virtual environment as the notebooks.
We provide the following scripts:

- [runRandomForest](scripts/runRandomForest.py): The script generates all of the ROC/PR plots that are presented in the manuscript and supplemental material.




## running the tool
A driver script is provided (``kce_tool.py``) as well as Jupyter notebooks that demonstrate the usage of the package.
``kce_tool.py`` has a number of commands that are used to implement the different functionalities. Run the script
with no arguments to see the commands. Run ``python kce_tool <command> -h`` to see the arguments for individual commands.

```
$ python kce_tool.py 
usage: kcet <command> [<args>]

The kcet commands are:
   pkpki        list of protein kinases and their inhibitors
   pkilist      get list of all protein kinase inhibitors

kinase cancer embedding tool

positional arguments:
  command     Subcommand to run

optional arguments:
  -h, --help  show this help message and exit
```

## Generate list of protein kinase inhibitors and corresponding protein kinases from DrugCentral
We create the file ``input/drug_kinase_links.tsv``which is obtained by applying the affinity(multiplicity) threshold 0.03
on data from DrugCentral. The file ``input/drug_kinase_links.tsv`` is a list of protein kinase inhibitors (PKI) that
have been used to treat cancer. Each PKI is shown together with its known major targets, act-value and a relevant
PubMed id (PMID).  To generate drug_kinase_links.tsv file, use the  following command:

``
python kce_tool.py pkpki [options]
``

``
--max_multiplicity: Limit on the number of PKs that are inhibited per PKI (default=5)
``

``
--outputfilename: The output file name, (default=input/drug_kinase_links.tsv)
``

## Generate list of protein kinase inhibitors
For the downstream analysis with the ``yactp`` tool, we need to have a file with the 
name of each PKI on one line and to remove duplicates. To generate this file, enter the following command
```
python kce_tool.py pkilist 
```

This script generates the file ``protein_kinase_inhibitors.txt``

## Use the yactp (yet another clinical trials parser) tool

Follow the instructions in the [yactp](https://github.com/monarch-initiative/yactp) repository
to generate the ``clinical_trials_by_phase.tsv`` file. We use the  ``protein_kinase_inhibitors.txt`` file
as the input for yactp. 
 

