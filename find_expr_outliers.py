#find gene expression outliers based on mean and std of mid expression values
#return file containing 0 if not outlier, 1 if yes

import pandas as pd
import sys

if len(sys.argv) != 7:
    print('python find_expr_outliers.py [gene expression file] '
    '[number of inds] [sd distance (e.g. 1,2,3)] [up/down/all (u/d/a)] '
    '[only loPoss genes? (y/n)] [out file]')
    #up/down/all refers to upper/higher expression outliers, down/lower expression
    #outliers or all outliers
    #loPoss refers to whether average population gene expression is too close to 0
    #so that there are no possible lower expression outliers
    sys.exit()

#command line arguments: expression file, output filename
exprpath = sys.argv[1]
n = float(sys.argv[2])
sdFactor = int(sys.argv[3])
oType = sys.argv[4] #outlier type: up/down/all
lPoss = sys.argv[5] #only genes where lower outlier is possible?
outpath = sys.argv[6]

#load gene expression data
gene_expr = pd.read_csv(exprpath, sep='\t')

#initialize new z-scores gene expression table
gene_expr_OL = pd.DataFrame(columns = gene_expr.columns)

for x in range(gene_expr.shape[0]):
    gene = gene_expr.iloc[x]
    geneSub = gene[6:]

    geneSub.sort_values(inplace=True)

    #sample of n individuals, use only those between 1st and 3rd quartile to
    #calculate mean and standard deviation
    st = int(n/4 + 0.5)
    en = int(n/4 * 3 + 0.5)
    geneInner = geneSub[st:en]

    mean = geneInner.mean()
    sd = geneInner.std()
    geneOL = []

    #only include genes where it's possible to have a lower outlier
    #so mean is greater than sdFactor * sd
    if (lPoss == 'y') & (mean < sdFactor*sd):
        continue

    #find outliers convert to 1, all others to 0
    for i in range(len(gene)):
        if i < 6:
            to_append = gene[i]
        elif oType == 'a':
            if abs(gene[i] - mean) > sdFactor*sd:
                to_append = 1
            else:
                to_append = 0
        elif oType == 'u':
            if (gene[i] - mean) > sdFactor*sd:
                to_append = 1
            else:
                to_append = 0
        elif oType == 'd':
            if (mean - gene[i]) > sdFactor*sd:
                to_append = 1
            else:
                to_append = 0
        geneOL.append(to_append)

    gene_series = pd.Series(geneOL, index = gene_expr_OL.columns)
    gene_expr_OL = gene_expr_OL.append(gene_series, ignore_index = True)

gene_expr_OL.to_csv(outpath, index=False, sep='\t')
