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
   pkpki        list of protein kinases and their inhibitors
   pkilist      get list of all protein kinase inhibitors

kinase cancer embedding tool

positional arguments:
  command     Subcommand to run

optional arguments:
  -h, --help  show this help message and exit
```

## Geenrate list of protein kinase inhibitors and corresponding protein kinases from DrugCentral
We create the file ``input/drug_kinase_links.tsv``which is obtained by applying the affinity(multiplicity) threshold 0.03
on data from DrugCentral. The file ``input/drug_kinase_links.tsv`` is a list of protein kinase inhibitors (PKI) that
have been used to treat cancer. Each PKI is shown together with its known major targets, act-value and a relevant
PubMed id (PMID). There is a limit on the number of PKs that are inhibited per PKI (default=5). To generate drug_kinase_links.tsv file, use the  following command:

``
 python kce_tool.py pkpki [options]
--max_multiplicity: Limit on the number of PKs that are inhibited per PKI (default=5)
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
as input for yactp.
 


# Running the notebooks

cd to the ``notebooks`` directory and enter
```
jupyter notebook
```
There are many notebooks to demonstrate and visualize data. There are 3 directories 
Predicting_links_all_phases, Predicting_links_phase_4 and 
Novel_Predictions which are for predicting  all phases using historical data,  phase 4 historical data and  novel predictions, respectively. 
To run the machine learning, first go to the notebook starting with ``DatSetGeneration`` and then go to the notebook starting with ``ExtractDifferenceVectors`` and follow directions there to download the required
embeddings and generate the difference vectors. Then go to the notebook starting with ``RandomForestPredictions``
to do RF learning. 
The directory data_exploration consists jupyter notebook files for visualization, showing
analogies and cosine similarities between words.

