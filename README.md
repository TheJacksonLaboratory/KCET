# kinaseCancerEmbeddingTool
Scripts and notebooks for the word-embedding kinase cancer prediction project:
In this repository, we provide instructions on how to obtain new kinase-cancer links.
See information for running the *Notebooks* at the bottom of this page.

## Set Up:
 Create a virtual environment in python 3 and activate the environment. 
```
virtualenv38
source p38/bin/activate
```

## running the tool
A driver script is provided (``kce_tool.py``) as well as Jupyter notebooks that demonstrate the usage of the package.
``kce_tool.py`` has a number of commands that are used to implement the different functionalities. Run the script
with no arguments to see the commands. Run ``python kce_tool <command> -h`` to see the arguments for individual commands.

```
$ python kce_tool.py 
usage: kcet <command> [<args>]

The kcet commands are:
   byphase      clinical trials by phase
   targeted     get list of targeted and untargeted kinases
   kinaselist   get list of all kinases
   merge        merge PKI/trial and PKI/kinase information
   ttp          write positive, negative, and prediction datasets

kinase cancer embedding tool

positional arguments:
  command     Subcommand to run

optional arguments:
  -h, --help  show this help message and exit
```




## Generate list of protein kinase inhibitors
The file ``input/drug_kinase_links.tsv`` is a hand-curated list of protein kinase inhibitors (PKI) that
have been used to treat cancer. Each PKI is shown together with its known major targets and a relevant
PubMed id (PMID). For the downstream analysis with the ``yactp`` tool, we need to have a file with the 
name of each PKI on one line and to remove duplicates. To generate this file, enter the following command

```
python kce_tool.py kinaselist
```
This script generates the file ``trainingset_protein_kinases.txt``

## Use the yactp (yet another clinical trials parser) tool

Follow the instructions in the [yatcp](https://github.com/monarch-initiative/yactp) repository
to generate the ``clinical_trials_by_phase.tsv`` file. We use the  ``trainingset_protein_kinases.txt`` file
as input for yactp.
 
## Collect protein kinase inhibitor trials by phase:
This script outputs two files; one represents the phase IV studies that have been published up to the
year of the -y argument. The other consists of all studies that
have been published up to the indicated year. If the ``-y`` argument is not used, then the
current year is taken.

```
python kce_tool.py byphase -c <path-to-clinical-trials-file> \
 [-y <start year>] \
 [--prefix]
```

The purpose of this script is severalfold. Our machine-learning experiment seeks to link
protein kinases (represented by gene symbols) to cancers that are treated by inhibiting the
kinase in question. The clinical trials data has links between protein kinase inhibitors 
and cancers. The file ``drug_kinase)_links.tsv`` has our hand-curated links between PKIs and the
kinases that are their main targets. The file ``prot_kinase.tsv`` further more links them
to the NCBI Gene ids for the kinases. 


* The clinical_trials_by_phase.tsv is required
* The ``-y`` argument is optional and defaults to the current year
* The ``--prefix`` argument is optional and defaults to KCET (kinase-cancer embedding tool). 

The output files are named based on the prefix. For example, if the default prefix is not changed and ``-y`` is set to 2018, the files will be called:

1. KCET_phaseIV_2018.tsv
2. KCET_all_phases_2018.tsv

Example:
```
python parse_ct_by_phase.py -c <path-to-clinical-trials-file>  -y 2018
```

## Find untargeted kinases

For training, we need to have a list of targeted kinases (positive examples) and untargeted kinases (negative examples).
In some cases, we will take only those kinases with phase 4 studies to be our positive examples.
Run  the script as follows:

```
python kce_tool.py targeted -c <path-to-clinical-trials-file>  [-y <year>]
```

This script outputs lists of targeted and untargeted kinases, as well as a list of kinasaes targeted by phase-4 studies.
For instance, if the ``-y`` argument is set to 2019, the following files are generated.


* targeted_kinases_2019.tsv 
* untargeted_kinases_2019.tsv
* targeted_kinases_phase_4_2019.tsv


## Generate positive/negative kinase-cancer links for training and testing

The links that we need for machine learning are between protein kinases and cancers. That is, if a protein kinase inhibitor (PKI)
has been successfully used to treat cancer X (phase 4 study started means post market approval for the PKI), and the PKI
inhibits kinases A and B, then we would generate the following positive dataset

```
A   X
B   X
```

We generate three files.

### positive

These files have links from phase 4 studies (i.e., the connection between the protein kinase inhibitor and the cancer was confirmed in other studies so that a post marketing, i.e., phase 4 study was started).

### negative

These files have ``factor`` times as many negative examples as we have positive examples. We do not use any example here that has
any clinical trial (phase 1, 2, 3, or 4), because even if phase 1 is not confirmatory, it does indicate there is some reason
to believe the PKI may be active against the cancer.

### prediction

This file contains all links except those in the positive set (the negatives are also included, since we believe that at least some of the negative links are actually positive).

To run the script, you do not need to previously run any of the other scripts (which are intended for visualization). 

```
python kce_tool.py ttp -c <path-to-clinical-trials-file> [-y <int>] [-f <int>] [-p <string>]
```

The ``-f`` argument is for the ``factor``, which is set to 10 by default. THe ``-p`` argument is the prefix, which is set to KCET by default.
Running the script for 2020 will produce these three files.

* KCET_positive_2020.tsv  
* KCET_negative_2020.tsv    
* KCET_prediction_2020.tsv


# Running the notebooks

cd to the ``noebooks`` directory and enter
```
jupyter notebook
```
There are many notebooks to demonstrate and visualize data. To run the machine learning, first go to the
notebook called ``ExtractDifferenceVectors.ipynb`` and follow directions there to download the required
embeddings and generate the difference vectors. Then go to the notebook ``RandomForestPredictions.ipynb``
to do RF learning. Finally, go to the notebook ``RankByDarkAndSLI.ipynb`` to do the final ranking.

