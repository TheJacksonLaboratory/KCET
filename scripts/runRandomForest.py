from kcet import KcetDatasetGenerator, KcetRandomForest

import pandas as pd
import numpy as np
import os
import sys
from sklearn.metrics import roc_auc_score, roc_curve, auc, precision_recall_curve, precision_score, recall_score, \
    average_precision_score, confusion_matrix
import matplotlib
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.abspath('..'))

plt.rc('axes', labelsize=18)
plt.rc('xtick', labelsize=16)  # fontsize of the tick labels
plt.rc('ytick', labelsize=16)  # fontsize of the tick labels

download_dir = '/home/peter/data/pubmed2vec'

ctfile = os.path.join(download_dir, "clinical_trials_by_phase.tsv")
embeddings2010 = os.path.join(download_dir, "embedding_SG_dim100_upto2010.npy")
words2010 = os.path.join(download_dir, "words_SG_upto2010.txt")
embeddings2014 = os.path.join(download_dir, "embedding_SG_dim100_upto2014.npy")
words2014 = os.path.join(download_dir, "words_SG_upto2014.txt")
if not os.path.isfile(ctfile):
    raise FileNotFoundError("Could not find clinical trials file at %s" % ctfile)
if not os.path.isfile(embeddings2010):
    raise FileNotFoundError("Could not find 2010 embeddings file at %s" % embeddings2010)
if not os.path.isfile(words2010):
    raise FileNotFoundError("Could not find 2010 words file at %s" % words2010)
if not os.path.isfile(embeddings2014):
    raise FileNotFoundError("Could not find 2014 embeddings file at %s" % embeddings2014)
if not os.path.isfile(words2014):
    raise FileNotFoundError("Could not find 2014 words file at %s" % words2014)


def year_label(begin_year: int, end_year: int):
    if begin_year == end_year:
        return "%d" % begin_year
    else:
        return "%d-%d" % (begin_year, end_year)


def plot_one_auc_curve(axis, y_test, yproba, begin_year: int, end_year: int, n_pos_test):
    fpr, tpr, thresholds_auc = roc_curve(y_test, yproba)
    auc_roc = roc_auc_score(y_test, yproba)
    yearl = year_label(begin_year, end_year)
    axis.plot(fpr, tpr, label='%s (%0.2f), n=%d' % (yearl, auc_roc, n_pos_test))
    return auc_roc


def plot_one_precision_recall_curve(axis, y_test, y_pred, begin_year: int, end_year: int, n_pos_test):
    """
    Plot a single precision recall curve
    axis: a matplotlib axis
    y_test: a numpy.ndarray with known classes
    y_pred: a numpy.ndarray with predictions
    midyear: an integer -- the start year for our predictions
    num_years_later: an integer -- how many years after midyear to go
    n_pos_test: number of true examples that are positive
    """
    precision, recall, thresholds = precision_recall_curve(y_test, y_pred)
    auc_recall_precision = average_precision_score(y_test, y_pred)
    yearl = year_label(begin_year, end_year)
    my_label = '%s, (%0.2f), n=%d' % (yearl, auc_recall_precision, n_pos_test)
    axis.plot(recall, precision, label=my_label)
    f1_scores = 2 * recall * precision / (recall + precision)
    # some nan values are encountered. The following replaces NaN by 0.0
    f1_scores = np.nan_to_num(f1_scores)
    best_threshold = thresholds[np.argmax(f1_scores)]
    best_f1 = np.max(f1_scores)
    precision_at_threshold = precision_score(y_test, y_pred > best_threshold)
    recall_at_threshold = recall_score(y_test, y_pred > best_threshold)
    return best_threshold, best_f1, precision_at_threshold, recall_at_threshold


def rrf(targetyear: int, test_years: list, outname: str, n_pk: int, phase4: bool):
    if targetyear == 2010:
        embeddings = embeddings2010
        words = words2010
    elif targetyear == 2014:
        embeddings = embeddings2010
        words = words2010
    else:
        raise ValueError("Invalid target year {}".format(targetyear))
    datagen = KcetDatasetGenerator(clinical_trials=ctfile, embeddings=embeddings, words=words, n_pk=n_pk)
    krf = KcetRandomForest(data_gen=datagen, target=targetyear, embedddingfile=embeddings2010, wordsfile=words2010)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    font = {'family': 'normal', 'size': 18}
    matplotlib.rc('font', **font)
    pr_data = []
    for (begin_y, end_y) in test_years:
        print(".", end='')
        y_pred, y_test, yproba, n_pos_test = krf.classify(begin_year=begin_y, end_year=end_y, phase4=phase4)
        auc_roc = plot_one_auc_curve(ax1, y_test, yproba, begin_y, end_y, n_pos_test)
        thresh, fscore, precision_at_threshold, recall_at_threshold = plot_one_precision_recall_curve(ax2, y_test,
                                                                                                      yproba, begin_y,
                                                                                                      end_y, n_pos_test)
        if phase4:
            phase = 'phase4'
        else:
            phase = 'allphases'
        pr_data.append({"target": targetyear,
                        "start": begin_y,
                        "end": end_y,
                        "n_pk": n_pk,
                        "phase": phase,
                        "AUROC": auc_roc,
                        "threshold": thresh,
                        "f-score": fscore,
                        "precision@threshold": precision_at_threshold,
                        "recall@threshold": recall_at_threshold,
                        "n_pos_test": n_pos_test})
    ax1.set_xlabel('1-Specificity');
    ax1.set_ylabel('Sensitivity');
    ax1.legend(loc="lower right");

    ax1.xaxis.set_tick_params(width=3)
    ax1.yaxis.set_tick_params(width=3)
    ax2.set_xlabel('Recall');
    ax2.set_ylabel('Precision');
    ax2.legend(loc="lower left");
    fig.savefig(outname, type='PDF')
    return pr_data


targetyear = 2010
prediction_years_by2 = [(2011, 2012), (2013, 2014), (2015, 2016), (2017, 2018), (2019, 2020)]
prediction_years_upto2020 = [(2011, 2011), (2011, 2012), (2011, 2014), (2011, 2016), (2011, 2018), (2011, 2020)]


def run_rrf(year: int):
    pr_data = []
    for npk in [1, 2, 5, 10]:
        outname = "m{}_{}_by_two_allphases.pdf".format(npk, year)
        dat = rrf(targetyear=year, test_years=prediction_years_by2, outname=outname, n_pk=npk, phase4=False)
        pr_data.append(dat)
        outname = "m{}_{}_by_two_phase4.pdf".format(npk, year)
        dat = rrf(targetyear=year, test_years=prediction_years_by2, outname=outname, n_pk=npk, phase4=True)
        pr_data.append(dat)
        outname = "m{}_{}_allyears_allphases.pdf".format(npk, year)
        dat = rrf(targetyear=year, test_years=prediction_years_upto2020, outname=outname, n_pk=npk, phase4=False)
        pr_data.append(dat)
        outname = "m{}_{}_allyears_phase4.pdf".format(npk, year)
        dat = rrf(targetyear=year, test_years=prediction_years_upto2020, outname=outname, n_pk=npk, phase4=True)
        pr_data.append(dat)
    return pr_data


prdat = run_rrf(2010)
df1 = pd.DataFrame.from_records(prdat)
df1.to_csv('pr2010.csv', index=False)
prdat = run_rrf(2014)
df1 = pd.DataFrame.from_records(prdat)
df1.to_csv('pr2014.csv', index=False)