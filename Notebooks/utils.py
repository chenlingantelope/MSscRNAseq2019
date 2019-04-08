import numpy as np
import matplotlib
save_path = '../CSF/Notebooks/'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import igraph as ig
import louvain
from sklearn.neighbors import kneighbors_graph
from random import sample
from numpy.random import permutation

def ClusterProp(clusters,batchid,cluster_description,tissue,states,celltypes,filename=None,clusterid=0):
    if clusterid==0:
        clusterid = np.unique(clusters[clusters>=0])
    nrows = (len(clusterid)//3)+1
    prop = []
    for i in np.unique(batchid):
        count = [np.sum(clusters[batchid==i]==j) for j in clusterid]
        prop.append(count/ np.sum(batchid==i))
    prop = np.asarray(prop)
    CSF = prop[(tissue=='CSF')]
    PBMC = prop[(tissue=='PBMC')]
    CSF_MS = prop[(tissue=='CSF')*(states=='MS')]
    CSF_Con = prop[(tissue=='CSF')*(states=='control')]
    PBMC_MS = prop[(tissue=='PBMC')*(states=='MS')]
    PBMC_Con = prop[(tissue=='PBMC')*(states=='control')]
    plt.figure(figsize=(12,3*nrows))
    for i,x in enumerate(clusterid):
        plt.subplot(nrows,3,i+1)
        boxes = plt.boxplot([CSF[:,i],PBMC[:,i],CSF_MS[:,i],CSF_Con[:,i],PBMC_MS[:,i],PBMC_Con[:,i]],patch_artist=True)
        colors=sns.color_palette("Paired",6)
        for patch, color in zip(boxes['boxes'], colors):
            patch.set_facecolor(color)
        plt.title(cluster_description[i] +' '+ celltypes[i] + " total %i cells"%sum(clusters==x))
        plt.xticks([1, 2, 3,4,5,6], ['CSF','PBMC','CSF MS', 'CSF Control', 'PBMC MS','PBMC Control'],rotation=45)
        plt.tight_layout()
    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)
    

def louvain_clusters(latent, k=10, rands=0, mutual=False):
    nn_matrix = kneighbors_graph(latent, k)
    rows, cols = nn_matrix.nonzero()
    if mutual == True:
        edges = [[row, col] if row < col else (col, row) for row, col in zip(rows, cols)]
        edges = np.asarray(edges)
        unique_edges, edges_count = np.unique(edges, return_counts=True, axis=0)
        edges = unique_edges[edges_count == 2]
    else:
        edges = [(row, col) for row, col in zip(rows, cols)]
    g = ig.Graph()
    g.add_vertices(latent.shape[0])
    g.add_edges(edges)
    louvain.set_rng_seed(rands)
    res = louvain.find_partition(g, louvain.ModularityVertexPartition)
    clusters = np.asarray(res.membership)
    return clusters


def matchedlist(norm_X, genelist, genenames):
    mean_exprs = np.asarray(norm_X[:,:].mean(axis=0)).ravel()
    plt.hist(np.log10(1+mean_exprs))
    matched=[]
    for x in genelist:
        if x in genenames:
            exprs = np.mean(norm_X[:,genenames==x])
            plt.axvline(np.log10(1+exprs),color='r')
            diff = np.abs(mean_exprs - exprs)
            idx = np.argsort(np.argsort(diff))
            idx = np.where((idx>1) * (idx<20))[0]
#             idx = sample(list(idx),1)[0]
            matched.append(genenames[idx])
            exprs = np.mean(norm_X[:,idx])
            plt.axvline(np.log10(1+exprs),color='b')
#     matched = np.asarray(matched)
#     print(genenames[matched])
#     return list(genenames[matched])
    return matched


def DotPlot(norm_X, genenames, genelist, clusterlabel, labelnames=None, filt=None,dotsize=500, filename=None,save_path=save_path+'figures/'):
    exprs = [norm_X[:,genenames==x] for x in genelist]
    exprs = np.asarray(exprs).squeeze()
    exprs = pd.DataFrame(exprs.T, columns=genelist)
    nz = deepcopy(exprs>0)
    exprs = (exprs-exprs.mean(axis=0))/exprs.std(axis=0)
    if labelnames is None:
        labelnames = list(np.unique(clusterlabel))
    if filt is None:
        filt = np.repeat(True,len(exprs))
    exprs = exprs.loc[filt]
    clusterlabel = clusterlabel[filt]
    nz = nz.loc[filt]
    mean_exprs = []
    mean_nz = []
    for x in np.unique(clusterlabel):
        avg = exprs.loc[clusterlabel==x].mean(axis=0)
        mean_exprs.append(avg)
        avg_nz = (nz.loc[clusterlabel==x]).mean(axis=0)
        mean_nz.append(avg_nz)
    mean_exprs = pd.concat(mean_exprs,axis=1)
    mean_exprs.columns = np.unique(clusterlabel)
    df = mean_exprs.T
#     df = (df-df.min())/(df.max()-df.min())
    df = np.log10(df)
    mean_exprs = df.T
    mean_exprs['genenames'] = list(mean_exprs.index)
    temp1 = mean_exprs.melt(id_vars=['genenames'], value_vars=np.unique(clusterlabel))
    mean_nz = pd.concat(mean_nz,axis=1)
    mean_nz.columns = np.unique(clusterlabel)
    mean_nz['genenames'] = mean_nz.index
    temp2 = mean_nz.melt(id_vars=['genenames'], value_vars=np.unique(clusterlabel))
    temp1['clusterid'] = [labelnames.index(x) for x in temp1['variable']]
    temp1['geneid'] = [genelist.index(x) for x in temp1['genenames']]
#     genes, temp1['geneid'] = np.unique(temp1['genenames'],return_inverse=True)
#     clustername, temp1['clusterid'] = np.unique(temp1['variable'],return_inverse=True)
    plt.figure(figsize=(8,4))
    plt.subplot(1,2,1)
    plt.scatter(x=temp1['geneid'],y=temp1['clusterid'],s=temp2['value']*dotsize*100, c=temp1['value'])
    plt.xticks(ticks=np.arange(len(genelist)), labels=genelist, rotation=90)
    plt.yticks(ticks=np.arange(len(labelnames)), labels=labelnames)
    plt.tight_layout()
    plt.colorbar()
    plt.subplot(1,2,2)
    for area in [0.1, 0.5, 1]:
        plt.scatter([], [], c='k', s=area*dotsize*100,
                    label=str(area))
        plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Proportion of Expression')
    plt.axis('off')
    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()
    return(mean_exprs,mean_nz)


def DotPlotCompare(norm_X, genenames, genelist, clusterlabel, group, filt,dotsize=500,labelnames=None,
                   filename=None,save_path=save_path+'figures/',
                  height=5,width=5):
    exprs = [norm_X[:,genenames==x] for x in genelist]
    exprs = np.asarray(exprs).squeeze()
    exprs = exprs / exprs.max(axis=1).reshape(len(genelist),1)
    exprs = pd.DataFrame(exprs.T, columns=genelist)
    exprs = exprs.loc[filt]
    clusterlabel = clusterlabel[filt]
    group = group[filt]
    mean_exprs = []
    clustercol = []
    groupcol = []
    colname = []
    if labelnames is None: 
        labelnames = np.unique(clusterlabel)
    for x in labelnames:
        for y in np.unique(group):
            if np.sum((clusterlabel==x)&(group==y))>0:
                avg=exprs.loc[(clusterlabel==x)&(group==y)].mean(axis=0)
                mean_exprs.append(np.log10(1+avg))
                colname.append(str(x)+'_'+str(y))
                clustercol.append(x)
                groupcol.append(y)
    mean_exprs = pd.concat(mean_exprs,axis=1)
    mean_exprs.columns = colname
    mean_exprs['genenames'] = mean_exprs.index
    temp = mean_exprs.melt(id_vars=['genenames'], value_vars=colname)
    temp['group'] = [x.split('_')[1] for x in temp['variable']]
    temp['cluster'] = [x.split('_')[0] for x in temp['variable']]
    groupname,temp['groupid'] = np.unique(temp['group'],return_inverse=True)
    genes = np.unique(temp['genenames'])
    temp['clusterid'] = [labelnames.index(x) for x in temp['cluster']]
    temp['geneid'] = [genelist.index(x) for x in temp['genenames']]
#     np.unique(temp['genenames'],return_inverse=True)
    plt.figure(figsize=(width*2,height))
    plt.subplot(1,2,1)
    plt.scatter(x=temp.loc[temp['groupid']==1]['geneid'],
                y=temp.loc[temp['groupid']==1]['clusterid'],
                s=temp.loc[temp['groupid']==1]['value']*dotsize*100,
                c='m',alpha=0.5,label=groupname[1])
    plt.scatter(x=temp.loc[temp['groupid']==0]['geneid'],
                y=temp.loc[temp['groupid']==0]['clusterid'],
                s=temp.loc[temp['groupid']==0]['value']*dotsize*100,
                c='c',alpha=0.5,label=groupname[0])
    plt.xticks(ticks=np.arange(len(genes)), labels=genelist, rotation=90)
    plt.yticks(ticks=np.arange(len(labelnames)), labels=labelnames)
    plt.legend(bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.subplot(1,2,2)
    for area in [0.001, 0.005, 0.01]:
        plt.scatter([], [], c='k', s=area*dotsize*100,
                    label=str(area))
        plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='log10 Expression')
    plt.axis('off')
    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()
    return temp


# def Heatmap(count,total,rankby,rownames,colnames,title,filename,width=9,height=8, save_path = save_path+'figures/'):
def Heatmap(count,total,rownames,colnames,title,filename,width=9,height=8, save_path = save_path+'figures/'):
    freq=[]
    nfreq = []
    for i in range(len(count[0,:])):
        f = (count[:,i]+1)/total
        freq.append(f)
        nf = np.mean(f)
        nfreq.append(np.log10(f/nf))

    freq = np.asarray(freq)
    nfreq = np.asarray(nfreq)
#     ranked = np.argsort(rankby)
    fig, ax = plt.subplots(figsize=(width,height))
    # We want to show all ticks...
#     plt.imshow(nfreq[ranked,:],aspect='auto',cmap='bwr')
    plt.imshow(nfreq,aspect='auto',cmap='bwr')
    plt.colorbar()
    ax.set_xticks(np.arange(len(rownames)))
    ax.set_yticks(np.arange(len(colnames)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(rownames,fontsize=20)
#     ax.set_yticklabels(colnames[ranked],fontsize=20)
    ax.set_yticklabels(colnames,fontsize=20)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    # Loop over data dimensions and create text annotations.
    for i in range(len(colnames)):
        for j in range(len(rownames)):
#             text = ax.text(j, i, "{:.0f}".format(count.T[ranked,:][i,j]),
            text = ax.text(j, i, "{:.0f}".format(count.T[i,j]),
                           ha="center", va="center",fontsize=15)
    ax.set_title(title,fontsize=30)
    fig.tight_layout()
    plt.savefig(save_path + filename, transparency=True)


