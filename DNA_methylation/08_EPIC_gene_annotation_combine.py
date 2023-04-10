import numpy as np
import statsmodels
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.stats.api as sms
import sys

EPIC_annot = []

for i in range(1,870001,10000): #for i in range(0,800000,10000):
    EPICnext = pd.read_csv("DNAm/data/annotation_split/EPIC_annotation_gene_" + str(i) + ".csv")
    EPIC_annot.append(EPICnext)


EPIC_annot = pd.concat(EPIC_annot)
EPIC_annot.to_csv("DNAm/data/EPIC_ensembl_gene_annotation.csv") 
