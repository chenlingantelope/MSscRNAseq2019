from scipy.stats import mannwhitneyu,wilcoxon
import numpy as np
from scipy.io import mmread
import pandas as pd

X = mmread('RFiles/all_data.mtx')
X = X.tocsr()

celllabels = np.load('Notebooks/meta/celllabels.npy')
isCSF = np.load('Notebooks/meta/isCSF.npy')
isMS = np.load('Notebooks/meta/isMS.npy')


logX = np.log10(1+X.todense())
scaling_factor = logX.mean(axis=1)
norm_X = logX - scaling_factor.reshape(len(scaling_factor), 1)

# def MannWhitneyUTest(norm_X, idx1, idx2):
#     res = []
#     for i in range(X.shape[1]):
#         x= np.asarray(X[idx1,i].todense()).ravel()
#         y= np.asarray(X[idx2,i].todense()).ravel()
#         if(len(np.unique(np.concatenate([x,y])))==1):
#             res.append([-1,-1])
#         else:
#             res.append(mannwhitneyu(x,y,alternative = 'two-sided'))
#     stat = np.asarray([x[0] for x in res])
#     pvalue = np.asarray([x[1] for x in res])
#     return(stat,pvalue)

def MannWhitneyUTest(X, idx1, idx2):
    res = []
    for i in range(X.shape[1]):
        x= np.asarray(X[idx1,i]).ravel()
        y= np.asarray(X[idx2,i]).ravel()
        if(len(np.unique(np.concatenate([x,y])))==1):
            res.append([-1,-1])
        else:
            res.append(mannwhitneyu(x,y,alternative = 'two-sided'))
    stat = np.asarray([x[0] for x in res])
    pvalue = np.asarray([x[1] for x in res])
    return(stat,pvalue)



celltypes = ['B1', 'B2', 'CD4', 'CD8a', 'CD8n', 'Gran', 'MegaK', 'Mono', 'NK1',
       'NK2', 'Tdg', 'Tregs', 'mDC1', 'mDC2', 'ncMono', 'pDC', 'plasma']


for i in celltypes:
    idx1 = (celllabels==i) & isCSF & isMS
    idx2 = (celllabels==i) & isCSF & (isMS == False)
    if (np.sum(idx1)>10) and (np.sum(idx2)>10):
        stat,pvalue = MannWhitneyUTest(norm_X, idx1, idx2)
        clusterid = np.repeat(i, len(stat))
        res = pd.DataFrame([clusterid,stat,pvalue],index=['clusterid','stat','pvalue']).T
        res.to_csv('DE/wilcoxon/MannWhitneyU.norm.MSinCSF.%s.csv'%i)


for i in celltypes:
    idx1 = (celllabels==i) & isCSF & (isMS == False)
    idx2 = (celllabels==i) & (isCSF == False) & (isMS == False)
    if (np.sum(idx1)>10) and (np.sum(idx2)>10):
        stat,pvalue = MannWhitneyUTest(norm_X, idx1, idx2)
        clusterid = np.repeat(i, len(stat))
        res = pd.DataFrame([clusterid,stat,pvalue],index=['clusterid','stat','pvalue']).T
        res.to_csv('DE/wilcoxon/MannWhitneyU.norm.tissue_control.%s.csv'%i)

for i in celltypes:
    idx1 = (celllabels==i) & (isCSF == False) & isMS
    idx2 = (celllabels==i) & (isCSF == False) & (isMS == False)
    if (np.sum(idx1)>10) and (np.sum(idx2)>10):
        stat,pvalue = MannWhitneyUTest(norm_X, idx1, idx2)
        clusterid = np.repeat(i, len(stat))
        res = pd.DataFrame([clusterid,stat,pvalue],index=['clusterid','stat','pvalue']).T
        res.to_csv('DE/wilcoxon/MannWhitneyU.norm.MSinPBMC.%s.csv'%i)

for i in celltypes:
    idx1 = (celllabels==i) & isMS
    idx2 = (celllabels==i) & (isMS == False)
    if (np.sum(idx1)>10) and (np.sum(idx2)>10):
        stat,pvalue = MannWhitneyUTest(norm_X, idx1, idx2)
        clusterid = np.repeat(i, len(stat))
        res = pd.DataFrame([clusterid,stat,pvalue],index=['clusterid','stat','pvalue']).T
        res.to_csv('DE/wilcoxon/MannWhitneyU.norm.MS.%s.csv'%i)

for i in celltypes:
    idx1 = (celllabels==i)
    idx2 = (celllabels!=i)
    if (np.sum(idx1)>10) and (np.sum(idx2)>10):
        stat,pvalue = MannWhitneyUTest(norm_X, idx1, idx2)
        clusterid = np.repeat(i, len(stat))
        res = pd.DataFrame([clusterid,stat,pvalue],index=['clusterid','stat','pvalue']).T
        res.to_csv('DE/wilcoxon/MannWhitneyU.norm.allclusters.%s.csv'%i)


norm_X = norm_X[celllabels=='CD4',:]
batchid = batchid[celllabels=='CD4']
isCSF = isCSF[celllabels=='CD4']
isMS = isMS[celllabels=='CD4']
celllabels = np.load('Notebooks/meta/CD4.clusters.npy')
celltypes = np.unique(celllabels)


for i in celltypes:
    idx1 = (celllabels==i) & isCSF & isMS
    idx2 = (celllabels==i) & isCSF & (isMS == False)
    if (np.sum(idx1)>10) and (np.sum(idx2)>10):
        stat,pvalue = MannWhitneyUTest(norm_X, idx1, idx2)
        clusterid = np.repeat(i, len(stat))
        res = pd.DataFrame([clusterid,stat,pvalue],index=['clusterid','stat','pvalue']).T
        res.to_csv('DE/wilcoxon/MannWhitneyU.CD4.MSinCSF.%s.csv'%i)


for i in celltypes:
    idx1 = (celllabels==i) & isCSF & (isMS == False)
    idx2 = (celllabels==i) & (isCSF == False) & (isMS == False)
    if (np.sum(idx1)>10) and (np.sum(idx2)>10):
        stat,pvalue = MannWhitneyUTest(norm_X, idx1, idx2)
        clusterid = np.repeat(i, len(stat))
        res = pd.DataFrame([clusterid,stat,pvalue],index=['clusterid','stat','pvalue']).T
        res.to_csv('DE/wilcoxon/MannWhitneyU.CD4.tissue_control.%s.csv'%i)

for i in celltypes:
    idx1 = (celllabels==i) & (isCSF == False) & isMS
    idx2 = (celllabels==i) & (isCSF == False) & (isMS == False)
    if (np.sum(idx1)>10) and (np.sum(idx2)>10):
        stat,pvalue = MannWhitneyUTest(norm_X, idx1, idx2)
        clusterid = np.repeat(i, len(stat))
        res = pd.DataFrame([clusterid,stat,pvalue],index=['clusterid','stat','pvalue']).T
        res.to_csv('DE/wilcoxon/MannWhitneyU.CD4.MSinPBMC.%s.csv'%i)

for i in celltypes:
    idx1 = (celllabels==i) & isMS
    idx2 = (celllabels==i) & (isMS == False)
    if (np.sum(idx1)>10) and (np.sum(idx2)>10):
        stat,pvalue = MannWhitneyUTest(norm_X, idx1, idx2)
        clusterid = np.repeat(i, len(stat))
        res = pd.DataFrame([clusterid,stat,pvalue],index=['clusterid','stat','pvalue']).T
        res.to_csv('DE/wilcoxon/MannWhitneyU.CD4.MS.%s.csv'%i)

for i in celltypes:
    idx1 = (celllabels==i)
    idx2 = (celllabels!=i)
    if (np.sum(idx1)>10) and (np.sum(idx2)>10):
        stat,pvalue = MannWhitneyUTest(norm_X, idx1, idx2)
        clusterid = np.repeat(i, len(stat))
        res = pd.DataFrame([clusterid,stat,pvalue],index=['clusterid','stat','pvalue']).T
        res.to_csv('DE/wilcoxon/MannWhitneyU.CD4.allclusters.%s.csv'%i)

