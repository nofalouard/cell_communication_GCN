import itertools
import random
import numpy as np
import pandas as pd
from itertools import compress
from itertools import combinations
from joblib import Parallel, delayed
from scipy import stats
import scanpy as sc
from anndata import AnnData
from matplotlib import pyplot as plt

def identify_interactions(adj_df, cell_annot):
    print("here1")
    inds_adj=adj_df.columns[0]
    adj_df.set_index(inds_adj, inplace=True)
    inds_cells=cell_annot.columns[0]
    cell_annot.set_index(inds_cells, inplace=True)
    print("here2")
    cells = np.asarray(cell_annot.index)
    labels = np.array(cell_annot['assign.level3_anno'], dtype = 'str')
    shuff = np.copy(labels)
    ctypes = list(set(labels)) 
    print("here3")
    label_dict = {ctypes[i]: list(compress(cells, labels== ctypes[i])) for i in range(len(ctypes))}
    cpairs = list(combinations(ctypes, 2)) + [(x,x) for x in ctypes]
    obs_means = {cpairs[i]: mean_pair_score(adj_df, label_dict,cpairs[i]) for i in range(len(cpairs))}
    #  dictionary stores null distribution (20 elements) of interaction scores between each pair of cell types
    null_dist = {cpairs[i]: calc_null_dist(shuff, cells, cpairs[i], adj_df) for i in range(len(cpairs))} 
    null_all =list(itertools.chain.from_iterable(null_dist.values()))
    pair_df = pd.DataFrame(index = ctypes, columns = ctypes)
    
    print("finished getting pair scores")
    '''
    for each cell type pair, calculate the p value of its observed mean score based on 
    its corresponding null distribution stored in the null_dist dictionary

    if observed mean > null distribution mean, then transform p value with -log10(pvalue)
    if observed mean < null distribution mean, then transform p value with log10(pvalue)

    if the observed mean is non-significant, i.e p value > 0.05, assign NA to its respective elements in the dataframe
    '''
    print("getting significance scores")
    for pair in cpairs:
        obs = obs_means[pair]
        null = null_dist[pair]
        dist = np.append(obs, np.copy(null))
        z_score = stats.zscore(dist)
        #print(z_score)
        if z_score[0] > 0:
            sig_scores = -np.log10(stats.norm.sf(abs(z_score)))
        else:
            sig_scores = np.log10(stats.norm.sf(abs(z_score)))
        if stats.norm.sf(abs(z_score))[0] >= 0.05:
            pair_df.loc[pair[0], pair[1]] = np.nan
            pair_df.loc[pair[1], pair[0]] = np.nan
        else:
            pair_df.loc[pair[0], pair[1]] = sig_scores[0] 
            pair_df.loc[pair[1], pair[0]] = sig_scores[0] 
    pair_df = pair_df.convert_dtypes()

    return pair_df

def interactions_by_clust(adj_df, cell_annot, clust_df):
    
    '''the commented code directly below is only needed if the dataframe arguments were read
       from csv files'''
    # print("here1")
    # inds_adj=adj_df.columns[0]
    # adj_df.set_index(inds_adj, inplace=True)
    # inds_cells=cell_annot.columns[0]
    # cell_annot.set_index(inds_cells, inplace=True)
    # clust_df = clust_df.drop(columns=clust_df.columns[0])
    clust1 = list(clust_df['Cluster_0'].dropna())
    clust2 = list(clust_df['Cluster_1'].dropna())
    # print("here2")
    # cells = np.asarray(cell_annot.index)
    cells = np.asarray(clust1+clust2)
    clust1_labels = np.array(cell_annot.loc[clust1, 'assign.level3_anno'], dtype = 'str')
    clust2_labels = np.array(cell_annot.loc[clust2, 'assign.level3_anno'], dtype = 'str')
    labels = np.array(cell_annot.loc[cells, 'assign.level3_anno'], dtype = 'str')
    shuff = np.copy(labels)
    shuff1 = np.copy(clust1_labels)
    shuff2 = np.copy(clust2_labels)
    ctypes1 = list(set(clust1_labels)) 
    ctypes2 = list(set(clust2_labels))
    ctypes = list(set(labels)) 
    label_dict1 = {ctypes1[i]: list(compress(clust1,clust1_labels== ctypes1[i])) for i in range(len(ctypes1))}
    label_dict2 = {ctypes2[i]: list(compress(clust2,clust2_labels== ctypes2[i])) for i in range(len(ctypes2))}
    cpairs1 = list(combinations(ctypes1, 2)) + [(x,x) for x in ctypes1]
    cpairs2 = list(combinations(ctypes2, 2)) + [(x,x) for x in ctypes2]
    cpairs = list(combinations(ctypes, 2)) + [(x,x) for x in ctypes]
    # print(cpairs)
    obs_means1 = {cpairs1[i]: mean_pair_score(adj_df, label_dict1,cpairs1[i]) for i in range(len(cpairs1))}
    obs_means2 = {cpairs2[i]: mean_pair_score(adj_df, label_dict2,cpairs2[i]) for i in range(len(cpairs2))}
    null_dist1 = {cpairs1[i]: calc_null_dist(shuff1, clust1, cpairs1[i], adj_df) for i in range(len(cpairs1))} 
    null_dist2 = {cpairs2[i]: calc_null_dist(shuff2, clust2, cpairs2[i], adj_df) for i in range(len(cpairs2))}
    null_dist = {cpairs[i]: calc_null_dist(shuff, cells, cpairs[i], adj_df) for i in range(len(cpairs))} 
    null_all =list(itertools.chain.from_iterable(null_dist.values()))
    null_clust1 =list(itertools.chain.from_iterable(null_dist1.values()))
    null_clust2 = list(itertools.chain.from_iterable(null_dist2.values()))
    pair_df1 = pd.DataFrame(index = ctypes1, columns = ctypes1)
    pair_df2 = pd.DataFrame(index = ctypes2, columns = ctypes2)
    print("starting pair1 analysis")
    for pair in cpairs1:
        obs = obs_means1[pair]
        null = null_dist1[pair]
        dist = np.append(obs, np.copy(null))
        z_score = stats.zscore(dist)
        #print(z_score)
        if z_score[0] > 0:
            sig_scores = -np.log10(stats.norm.sf(abs(z_score)))
        else:
            sig_scores = np.log10(stats.norm.sf(abs(z_score)))
        if stats.norm.sf(abs(z_score))[0] >= 0.05:
            pair_df1.loc[pair[0], pair[1]] = np.nan
            pair_df1.loc[pair[1], pair[0]] = np.nan
        else:
            pair_df1.loc[pair[0], pair[1]] = sig_scores[0] 
            pair_df1.loc[pair[1], pair[0]] = sig_scores[0] 
    pair_df1 = pair_df1.convert_dtypes()
    print("starting pair2 analysis")
    for pair in cpairs2:
        obs = obs_means2[pair]
        null = null_dist2[pair]
        dist = np.append(obs, np.copy(null))
        z_score = stats.zscore(dist)
        #print(z_score)
        if z_score[0] > 0:
            sig_scores = -np.log10(stats.norm.sf(abs(z_score)))
        else:
            sig_scores = np.log10(stats.norm.sf(abs(z_score)))
        if stats.norm.sf(abs(z_score))[0] >= 0.05:
            pair_df2.loc[pair[0], pair[1]] = np.nan
            pair_df2.loc[pair[1], pair[0]] = np.nan
        else:
            pair_df2.loc[pair[0], pair[1]] = sig_scores[0] 
            pair_df2.loc[pair[1], pair[0]] = sig_scores[0] 
    pair_df2 = pair_df2.convert_dtypes()

    return pair_df1, pair_df2
    
def ctype_ratios(cell_annot, clust_df):
    '''the commented code directly below is only needed if the dataframe arguments were read
       from csv files'''
    # inds_cells=cell_annot.columns[0]
    # cell_annot.set_index(inds_cells, inplace=True)
    # clust_df = clust_df.drop(columns=clust_df.columns[0])
    clust1 = list(clust_df['Cluster_0'].dropna())
    clust2 = list(clust_df['Cluster_1'].dropna())
    all_cells = clust1+clust2
    clust1_labels = np.array(cell_annot.loc[clust1, 'assign.level3_anno'], dtype = 'str')
    clust2_labels = np.array(cell_annot.loc[clust2, 'assign.level3_anno'], dtype = 'str')
    clust_labels = np.array(cell_annot.loc[all_cells,'assign.level3_anno'], dtype = 'str')
    ctypes1 = list(set(clust1_labels)) 
    ctypes2 = list(set(clust2_labels))
    ctypes = list(set(clust_labels))
    label_dict1 = {ctypes1[i]: list(compress(clust1,clust1_labels== ctypes1[i])) for i in range(len(ctypes1))}
    label_dict2 = {ctypes2[i]: list(compress(clust2,clust2_labels== ctypes2[i])) for i in range(len(ctypes2))}
    label_dict = {ctypes[i]: list(compress(all_cells,clust_labels== ctypes[i])) for i in range(len(ctypes))}
    ctype_count = {ctypes[i]: len(label_dict[ctypes[i]]) for i in range(len(ctypes))}
    clust1_ratios = {ctypes1[i]: len(label_dict1[ctypes1[i]])/ctype_count[ctypes1[i]]*100 for i in range(len(ctypes1))}
    clust2_ratios = {ctypes2[i]: len(label_dict2[ctypes2[i]])/ctype_count[ctypes2[i]]*100 for i in range(len(ctypes2))}
    ratio_compare = [clust1_ratios, clust2_ratios]
    # for c in ctypes1:
    #     print(clust1_ratios[c] + clust2_ratios[c])
    # # print(list(label_dict1.values()))
    # print(len(clust1), len(clust2), len(all_cells))
    # sum_clust1 = 0
    # sum_clust2 = 0
    # for lst in label_dict2.values():
    #     sum_clust2 += len(lst)
    # for lst in label_dict1.values():
    #     sum_clust1 += len(lst)
    # print(sum_clust1, sum_clust2, len(all_cells))
    # print(len(clust1_labels), len(clust_labels))
    # df_lst = [pd.DataFrame.from_dict(ctype_count, orient = 'index', columns = ['Cell_Type_Counts'])]
    df_lst = []
    for i, ratios in enumerate(ratio_compare):
        df = pd.DataFrame.from_dict(ratios, orient = 'index', columns = ['Cluster_' + str(i+1)])
        df_lst.append(df)
    merged_dfs = pd.concat(df_lst, axis = 1)

    print("merged dfs")

    fig, axs = plt.subplots(1, 1,figsize =(10, 7),tight_layout = True)
    plt.xlabel("Cell Type")
    plt.ylabel("Percentage")
    plt.title('Cell Type Percentages in Each Cluster')
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    # plt.xticks(fontsize=16)
    # plt.xticks(np.arange(0, 110, step=10))
    # complement_ratios = [100 - x for x in list(clust1_ratios.values())]
  
    axs.bar(list(clust1_ratios.keys()), [100] * len(list(clust1_ratios.keys())), label = 'Cluster 1')
    # axs.bar(list(clust1_ratios.keys()), list(clust1_ratios.values()), label = 'Cluster 1')
    axs.bar(list(clust2_ratios.keys()), list(clust2_ratios.values()), label = 'Cluster 2')
    axs.legend()
    return plt
    # return merged_dfs
'''
mean_pair_score:
Given df, a weighted adjacency matrix, calculates the mean value of edges lying between nodes (cells) that belong to
the two types specified in the "pair" argument
'''
def mean_pair_score(df, label_dict, pair):
    c1c2_mean = np.asarray(df.loc[label_dict[pair[0]], label_dict[pair[1]]]).mean()
    c2c1_mean = np.asarray(df.loc[label_dict[pair[1]], label_dict[pair[0]]]).mean()
    pair_mean = (c1c2_mean + c2c1_mean)/2
    return pair_mean

'''
calc_null_dist:
Given the list of cell type labels for each cell, randomly permutes the labels 1000x
to generate a null distribution of average interaction score for each pair of cell types

Returns the null distribution as a list of 1000 elements
'''
def calc_null_dist(shuff_labels, cell_names, pair_select, df):
    null_lst = []
#     print("in func")
    for i in range(20):
        random.shuffle(shuff_labels)
        shuff_dict = {pair_select[j]: list(compress(cell_names,shuff_labels== pair_select[j])) for j in range(2)}
        score = mean_pair_score(df, shuff_dict, pair_select)
        null_lst.append(score)
    return null_lst